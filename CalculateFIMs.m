function [Fs,Gs,ms,cons,Ts] = CalculateFIMs(cluster, modelgenfun, mfile, njobs, Ts, submitOpts, maxWait, useObj, useConstraints, varargin)

% Determine which parameter sets will be calculated by which jobs
np = size(Ts,2);
np_per_job = repmat(floor(np/njobs), njobs, 1);
np_left = np-sum(np_per_job);
np_per_job(1:np_left) = np_per_job(1:np_left)+1;
pi_per_job = mat2cell((1:np)', np_per_job, 1);

rankfun = @(F,G)[];

rngseed = 0;
nparamsetsperjob = []; % Doesn't matter what this is set to when supplying parameter sets

% Set defaults for cluster submission options
submitOpts = defaultSubmitOpts(submitOpts, 'CalculateFIMs');
nworkersperjob = submitOpts.proc;

inputs = cell(njobs,1);
for i = 1:njobs
    samplerfun_i = @(nparams,nsamps,lb,ub)Ts(:,pi_per_job{i});
    inputs{i} = [{nworkersperjob, @CalculateFIMs_job, modelgenfun, mfile, rankfun, nparamsetsperjob, samplerfun_i, useObj, useConstraints, rngseed} varargin(:)'];
end

submitOpts.name = 'CalculateFIMs';
if ~isfield(submitOpts, 'additionalFun')
    submitOpts.additionalFun = {};
end
submitOpts.additionalFun = [submitOpts.additionalFun(:); {mfile;'CalculateFIMs_job';func2str(modelgenfun)}];

nout = 5;
results = submitClusterSync(cluster, @InitializeParallel, inputs, nout, submitOpts, 60, maxWait);
[ranks,ms,cons,Gs,Fs] = deal(cell(njobs,1));
jobiscomplete = ~cellfun(@isempty, results);
jobcomplete_i = find(jobiscomplete);
for i = jobcomplete_i(:)'
    [ranks{i},ms{i},cons{i},Gs{i},Fs{i}] = deal(results{i}{:});
end
ranks = vertcat(ranks{:});
ms = vertcat(ms{:});
cons = [cons{:}];
Gs = vertcat(Gs{:});
Fs = cat(3, Fs{:});
Tskept = vertcat(pi_per_job{jobiscomplete});
Ts = Ts(:,Tskept);
        
end