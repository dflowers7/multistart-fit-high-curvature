function direc = DetermineGoodG_submit(cluster, njobs, additionalFun, jobnamesuffix, rngfactor, modelgenfun, varargin)

% Process inputs
p = inputParser();
p.addRequired('cluster');
p.addRequired('njobs');
p.addRequired('additionalFun');
p.addRequired('jobnamesuffix');
p.addRequired('rngfactor');
p.addRequired('modelgenfun');
p.addOptional('GRelTol', []);
p.addOptional('fitfun', []);
p.addOptional('useObj', []);
p.addOptional('useConstraints',[]);
p.addOptional('initopts',[]);
p.parse(cluster, njobs, additionalFun, jobnamesuffix, rngfactor, modelgenfun, varargin{:});
res = p.Results;
cluster = res.cluster;
njobs = res.njobs;
additionalFun = res.additionalFun;
jobnamesuffix = res.jobnamesuffix;
rngfactor = res.rngfactor;
modelgenfun = res.modelgenfun;
GRelTol = res.GRelTol;
fitfun = res.fitfun;
useObj = res.useObj;
useConstraints = res.useConstraints;
initopts = res.initopts;

workersperjob = 4;

inputs = cell(njobs,1);
for i = 1:njobs
    rngseed = i*rngfactor;
    inputs{i} = {workersperjob, @DetermineGoodG_job, modelgenfun, rngseed, GRelTol, fitfun, useObj, useConstraints, initopts};
end

nout = 5;
submitOpts.name = ['DetermineGoodG_' jobnamesuffix];
submitOpts.proc = 4;
submitOpts.time = 24;
submitOpts.jvm = true;
customAdditionalFun = additionalFun;
submitOpts.additionalFun = [{'DetermineGoodG_job';func2str(modelgenfun)}; customAdditionalFun(:)];

direc = submitCluster(cluster, @InitializeParallel, inputs, nout, submitOpts);

end