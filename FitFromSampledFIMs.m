function direc = FitFromSampledFIMs(cluster, FIMdirecs, Gstop, maxNoAttempts, modelgenfun_arginnames, keepms, keepcons, submitOpts, njobsperconfig, varargin)

if ischar(FIMdirecs)
    FIMdirecs = cellstr(FIMdirecs);
end

if iscellstr(FIMdirecs)
    % Create input argument tables in each directory to speed up loading of
    % FIMs
    inputargtables = cell(numel(FIMdirecs),1);
    JobFileOutputs = cell(numel(FIMdirecs),1);
    for di = 1:numel(FIMdirecs)
        inputargtables{di} = CreateInputArgTable(FIMdirecs{di});
        JobFileOutputs{di} = cell(6,1);
        [JobFileOutputs{di}{:}] = JobFiles(FIMdirecs{di});
    end
end

% Process modelgenfun input arguments
assert(all(cellfun(@iscell, varargin)), 'Each extra input argument to this function should be a cell array enumerating the possible inputs for each of modelgenfun''s input arguments.')
nargin_modelgenfun = numel(varargin);
npossibleargs_modelgenfun = cellfun(@numel, varargin);
nconfigs = prod(npossibleargs_modelgenfun);
njobs = nconfigs*njobsperconfig;

% Set defaults for cluster submission options
if isempty(submitOpts)
    submitOpts = struct;
end
defaultSubmitOpts.name = 'FitFromSampledFIMs';
defaultSubmitOpts.proc = 4;
defaultSubmitOpts.time = 24;
defaultSubmitOpts.jvm = true;
submitOpts = mergestruct(defaultSubmitOpts, submitOpts);

nworkersperjob = submitOpts.proc;

% Split the number of attempted fits up among the jobs per configuration
np_per_job = repmat(floor(maxNoAttempts/njobsperconfig), njobsperconfig, 1);
np_leftover = maxNoAttempts - sum(np_per_job);
np_per_job(1:np_leftover) = np_per_job(1:np_leftover) + 1;
pi_per_job = mat2cell(1:maxNoAttempts, 1, np_per_job);

% For each configuration of modelgenfun input args, find the relevant FIM
% information
inputs = cell(njobs,1);
modelgenfuns = cell(njobs,1);
mfiles = cell(njobs,1);
lastInputArgFilter = 'unset';
lasttab_ji = 'unset';
for ji = 1:njobs
    
    % Get the modelgenfun input arguments for this job
    whichargs = cell(nargin_modelgenfun,1);
    [pset_i,whichargs{:}] = ind2sub([njobsperconfig npossibleargs_modelgenfun], ji);
    varargin_i = struct;
    for j = 1:nargin_modelgenfun
        varargin_i.(modelgenfun_arginnames{j}) = varargin{j}{whichargs{j}};
    end
    % Create an input argument filter for these input arguments
    inputArgFilter = varargin_i;
    
    % Use the input argument filter to load only the FIMs used for this set
    % of input arguments
    if ~isequal(inputArgFilter, lastInputArgFilter)
        tab_ji = CollectSampledFIMs(FIMdirecs, modelgenfun_arginnames, keepms, keepcons, inputArgFilter, inputargtables, JobFileOutputs);
    else
        tab_ji = lasttab_ji;
    end
    lastInputArgFilter = inputArgFilter;
    lasttab_ji = tab_ji;
    
    % Sort the FIM ranks descending
    [rankssort,rankssort_i] = sort(vertcat(tab_ji.ranks), 1, 'descend');
    p_i = rankssort_i(pi_per_job{pset_i});
    
    % Collect parameters
    if keepms
        ms_ji = vertcat(tab_ji.ms);
        ks_ji = [ms_ji.k];
    else
        ks_ji = [];
    end
    if keepcons
        cons_ji = [tab_ji.cons];
        [ss_ji,qs_ji,hs_ji] = deal(cell(size(cons_ji,1),1));
        for coni = 1:size(cons_ji,1)
            ss_ji{coni} = [cons_ji(coni,:).s];
            qs_ji{coni} = [cons_ji(coni,:).q];
            hs_ji{coni} = [cons_ji(coni,:).h];
        end
    else
        ss_ji = [];
        qs_ji = [];
        hs_ji = [];
    end
    
    % Collect varargin
    varargin_ji = struct2cell(varargin_i);
    
    % Collect other inputs
    modelgenfuns{ji} = tab_ji(1).modelgenfun;
    mfiles{ji} = tab_ji(1).mfile;
    
    % Sort parameters by descending ranks, only keeping the top
    % maxNoAttempts results
    if keepms
        ks_ji = ks_ji(:,p_i);
    end
    if keepcons
        for coni = 1:size(cons_ji,1)
            ss_ji{coni} = ss_ji{coni}(:,p_i);
            qs_ji{coni} = qs_ji{coni}(:,p_i);
            hs_ji{coni} = hs_ji{coni}(:,p_i);
        end
    end
    
    inputs{ji} = [{nworkersperjob, @FitFromSampledFIMs_job, modelgenfuns{ji}, mfiles{ji}, ks_ji, ss_ji, qs_ji, hs_ji, Gstop} varargin_ji(:)'];
end

if ~isfield(submitOpts, 'additionalFun')
    submitOpts.additionalFun = {};
end
mfiles = unique(mfiles);
modelgenfuns = unique(cellfun(@func2str, modelgenfuns, 'UniformOutput', false));
submitOpts.additionalFun = [submitOpts.additionalFun(:); mfiles; modelgenfuns; {'FitFromSampledFIMs'}];
nout = 4;
direc = submitCluster(cluster, 'InitializeParallel', inputs, nout, submitOpts);

end

function [mbest, conbest, Gbest, Dbest] = FitFromSampledFIMs_job(modelgenfun, mfile, ks, ss, qs, hs, Gstop, varargin)

[m,con,obj,opts] = modelgenfun(mfile, varargin{:});

% Determine the number of fits to be run
p = {ks, ss, qs, hs};
pisempty = cellfun(@isempty, p);
if ~pisempty(1)
    nFits = size(ks,2);
else
    nFits_i = find(~pisempty, 1);
    nFits = size(p{nFits_i}{1}, 2);
end

ncon = numel(con);
Gbest = inf;
for i = 1:nFits
    
    if ~isempty(ks)
        m_i = m.Update(ks(:,i));
    else
        m_i = m;
    end
    
    con_i = con;
    for j = ncon:-1:1
        if isempty(ss)
            s_j = con_i(j).s;
        else
            s_j = ss{j}(:,i);
        end
        if isempty(qs)
            q_j = con_i(j).q;
        else
            q_j = qs{j}(:,i);
        end
        if isempty(hs)
            h_j = con_i(j).h;
        else
            h_j = hs{j}(:,i);
        end
        con_i(j) = con_i(j).Update(s_j, q_j, h_j);
    end
    
    [m_i,con_i,G_i,D_i] = FitObjective(m_i, con_i, obj, opts);
    
    if G_i < Gbest
        mbest = m_i;
        conbest = con_i;
        Gbest = G_i;
        Dbest = D_i;
    end
    
    if Gbest <= Gstop
        break
    end
end

end