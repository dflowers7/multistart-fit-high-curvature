function direc = FitAllSampledFIMs(cluster, FIMdirecs, njobsperconfig, modelgenfun_arginnames, keepms, keepcons, submitOpts, varargin)
% direc = FitAllSampledFIMs(cluster, FIMdirecs, njobsperconfig, modelgenfun_arginnames, keepms, keepcons, submitOpts, varargin)

FIMdirecs = cellstr(FIMdirecs);

% Create input argument tables for each directory to speed up loading of
% FIMs
inputargtables = cell(numel(FIMdirecs),1);
JobFileOutputs = cell(numel(FIMdirecs),1);
for di = 1:numel(FIMdirecs)
    inputargtables{di} = CreateInputArgTable(FIMdirecs{di});
    JobFileOutputs{di} = cell(6,1);
    [JobFileOutputs{di}{:}] = JobFiles(FIMdirecs{di});
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
defaultSubmitOpts.name = 'FitAllSampledFIMs';
defaultSubmitOpts.proc = 4;
defaultSubmitOpts.time = 24;
defaultSubmitOpts.jvm = true;
submitOpts = mergestruct(defaultSubmitOpts, submitOpts);

nworkersperjob = submitOpts.proc;

% For each configuration of modelgenfun input args, find the relevant FIM
% information
inputs = cell(njobs,1);
modelgenfuns = cell(njobs,1);
mfiles = cell(njobs,1);
lastInputArgFilter = 'unset';
lasttab_config = 'unset';
for ji = 1:njobs
    
    % Get the modelgenfun input arguments for this job
    whichargs = cell(nargin_modelgenfun,1);
    [jperconfig_i,whichargs{:}] = ind2sub([njobsperconfig npossibleargs_modelgenfun], ji);
    varargin_i = struct;
    for j = 1:nargin_modelgenfun
        varargin_i.(modelgenfun_arginnames{j}) = varargin{j}{whichargs{j}};
    end
    % Create an input argument filter for these input arguments
    inputArgFilter = varargin_i;
    
    % Use the input argument filter to load only the FIMs used for this set
    % of input arguments
    if ~isequal(inputArgFilter, lastInputArgFilter)
        tab_config = CollectSampledFIMs(FIMdirecs, modelgenfun_arginnames, keepms, keepcons, inputArgFilter, inputargtables, JobFileOutputs);
        tab_config = ExpandCollectedFIMs(tab_config);
    else
        tab_config = lasttab_config;
    end
    lastInputArgFilter = inputArgFilter;
    lasttab_config = tab_config;
       
    % If it is the first job for a set of modelgenfun input arguments,
    % split the parameter sets into their separate jobs
    if jperconfig_i == 1
        basenparamsetsperjob = floor(numel(tab_config)/njobsperconfig);
        nparamsetsperjob = repmat(basenparamsetsperjob, njobsperconfig, 1);
        remsets = numel(tab_config)-sum(nparamsetsperjob);
        nparamsetsperjob(1:remsets) = nparamsetsperjob(1:remsets)+1;
    end
    tab_ji = mat2cell(tab_config, nparamsetsperjob, 1);
    tab_ji = tab_ji{jperconfig_i};
    
    % Sort the FIM ranks descending
%     [rankssort,rankssort_i] = sort(vertcat(tab_ji.ranks), 1, 'descend');
%     topn_i = rankssort_i(1:maxNoAttempts);
    
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
    
    % Sort parameters by descending ranks, only keeping the top
    % maxNoAttempts results
%     if keepms
%         ks_ji = ks_ji(:,topn_i);
%     end
%     if keepcons
%         for coni = 1:size(cons_ji,1)
%             ss_ji{coni} = ss_ji{coni}(:,topn_i);
%             qs_ji{coni} = qs_ji{coni}(:,topn_i);
%             hs_ji{coni} = hs_ji{coni}(:,topn_i);
%         end
%     end
    
    % Collect varargin
    varargin_ji = struct2cell(varargin_i);
    
    % Collect other inputs
    modelgenfuns{ji} = tab_ji(1).modelgenfun;
    mfiles{ji} = tab_ji(1).mfile;
    
    inputs{ji} = [{nworkersperjob, @FitAllSampledFIMs_job, modelgenfuns{ji}, mfiles{ji}, ks_ji, ss_ji, qs_ji, hs_ji} varargin_ji(:)'];
end

if ~isfield(submitOpts, 'additionalFun')
    submitOpts.additionalFun = {};
end
mfiles = unique(mfiles);
modelgenfuns = unique(cellfun(@func2str, modelgenfuns, 'UniformOutput', false));
submitOpts.additionalFun = [submitOpts.additionalFun(:); mfiles; modelgenfuns; {'FitAllSampledFIMs'}];
nout = 4;
direc = submitCluster(cluster, 'InitializeParallel', inputs, nout, submitOpts);

end

function [ms, cons, Gs, Ds] = FitAllSampledFIMs_job(modelgenfun, mfile, ks, ss, qs, hs, varargin)

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
ms = struct('Name', cell(nFits,1), 'Parameters', cell(nFits,1));
cons = struct('Name', cell(ncon,nFits), 's', cell(ncon,nFits), 'q', cell(ncon,nFits), 'h', cell(ncon,nFits));
Gs = nan(nFits,1);
Ds = [];
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
    
    if isempty(Ds)
        % Once we know the number of fit parameters, initialize the
        % gradient store
        Ds = nan(numel(D_i), nFits);
    end
    
    ms(i) = struct('Name', m_i.Name, 'Parameters', m_i.Parameters);
    cons(:,i) = struct('Name', {con_i.Name}.', 's', {con_i.s}.', 'q', {con_i.q}.', 'h', {con_i.h}.');
    Gs(i) = G_i;
    Ds(:,i) = D_i;
end

end