function direc = SampleFIMs(cluster, modelgenfun, mfile, njobsperconfig, nFIMsperjob, sampleOpts, submitOpts, rngseed_start, rngseed_interval, varargin)
% direc = SampleFIMs(modelgenfun, mfile, nFIMs, nworkersperjob, sampleOpts, submitOpts, varargin)
% Submits jobs that sample parameters uniformly between the bounds and
% calculates their FIMs.
% Input arguments:
%   cluster
%   modelgenfun
%   mfile
%   njobsperconfig
%   nFIMSperjob
%   sampleOpts
%       .IncludeFIMsInOutput [false]
%       .UseLogSampling [true]
%       .UseConstraintObjs [ {true;true;...;true} ]
%       .UseObjs    [ all true ]
%   submitOpts
%       See submitCluster.m.
%   varargin
%       Possible values for other input arguments to modelgenfun go here.
%       Each should be a cell array of possible input arguments.

% Set defaults for FIM sampling options
if isempty(sampleOpts)
    sampleOpts = struct;
end
defaultSampleOpts.IncludeFIMsInOutput = false;
defaultSampleOpts.UseLogSampling = true;
defaultSampleOpts.UseConstraintObjs = [];
defaultSampleOpts.UseObjs = [];
defaultSampleOpts.RankFunction = @(F,G)calculateMetric(F,G,'EigenvalueThresholdCount',1e-3);
sampleOpts = mergestruct(defaultSampleOpts, sampleOpts);
includeFIMSinoutput = sampleOpts.IncludeFIMsInOutput;
uselogsampling = sampleOpts.UseLogSampling;
useConstraints = sampleOpts.UseConstraintObjs;
useObj = sampleOpts.UseObjs;
rankfun = sampleOpts.RankFunction;

% Set defaults for cluster submission options
if isempty(submitOpts)
    submitOpts = struct;
end
submitOpts = defaultSubmitOpts(submitOpts, 'SampleFIMs');

nworkersperjob = submitOpts.proc;

samplerfun = @(nT,nsamples,lowerbounds,upperbounds)uniformlySampleParameters(nT,nsamples,lowerbounds,upperbounds,uselogsampling);

% varargin should be a cell array of cell arrays
assert(all(cellfun(@iscell, varargin)), 'Each extra input argument to this function should be a cell array enumerating the possible inputs for each of modelgenfun''s input arguments.')
nargin_modelgenfun = numel(varargin);
npossibleargs_modelgenfun = cellfun(@numel, varargin);
njobs = prod(npossibleargs_modelgenfun)*njobsperconfig;

inputs = cell(njobs,1);
rngseed_i = rngseed_start;
for i = 1:njobs
    whichargs = cell(nargin_modelgenfun,1);
    [~,whichargs{:}] = ind2sub([njobsperconfig npossibleargs_modelgenfun], i);
    varargin_i = cell(nargin_modelgenfun,1);
    for j = 1:nargin_modelgenfun
        varargin_i{j} = varargin{j}{whichargs{j}};
    end
    inputs{i} = [{nworkersperjob, @CalculateFIMs_job, modelgenfun, mfile, rankfun, nFIMsperjob, samplerfun, useObj, useConstraints, rngseed_i}, varargin_i(:)'];
    rngseed_i = rngseed_i + rngseed_interval;
end

nout = 4 + includeFIMSinoutput;

if ~isfield(submitOpts, 'additionalFun')
    submitOpts.additionalFun = {};
end
submitOpts.additionalFun = [submitOpts.additionalFun(:); {mfile;'SampleFIMs';func2str(modelgenfun)}];
direc = submitCluster(cluster, 'InitializeParallel', inputs, nout, submitOpts);

end