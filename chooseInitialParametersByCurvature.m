function [k,Fs,p] = chooseInitialParametersByCurvature(m, con, obj, fitopts, initopts, uselogsampling)
% initopts
%   .UseHessian
%   .Metric
%       'SumOfConditionNumbers'
%       'EigenvalueThresholdCount'
%       'ConditionNumberThresholdCount'
%       'ObjectiveScaledEigenvalueThresholdCount'
%       'SumOfObjectiveScaledEigenvalues'
%       'Random'
%   .Threshold
%   .NumberOfSamples
% uselogsampling [ logical scalar or nT-length vector {true} ]
%   If a scalar, all parameters will be sampled uniformly on a log
%   basis if true and uniformly on a linear basis if false. If a vector,
%   one can indicate for each individual fit parameter whether it should be
%   sampled uniformly on a log (true) or linear (false) basis.

if nargin < 6
    uselogsampling = true;
end

nParamSets = 1; % This function is currently used to select one parameter set only, but I'm leaving this here in case I want to add functionality

% Calculate the number of fit parameters
nTk = sum(fitopts.UseParams);
nTs = sum(fitopts.UseSeeds(:));
nTq = sum(cellfun(@sum, fitopts.UseInputControls));
nTh = sum(cellfun(@sum, fitopts.UseDoseControls));

nT = nTk + nTs + nTq + nTh;
pLo = fitopts.LowerBound;
pHi = fitopts.UpperBound;

% Standardize uselogbasis as a vector and check its type and length
if isscalar(uselogsampling)
    uselogsampling = repmat(uselogsampling, nT, 1);
end
assert(islogical(uselogsampling), 'uselogsampling should be a logical scalar or vector.')
assert(numel(uselogsampling) == nT, ...
    'uselogsampling had %g elements but should have %g elements, the number of parameters being fit.', ...
    numel(uselogsampling), nT)

if isfield(initopts, 'NumberOfSamples')
    nParamTrys = initopts.NumberOfSamples;
else
    nParamTrys = [];
end
if isempty(nParamTrys)
    nParamTrys = 100*nParamSets;
end
p = uniformlySampleParameters(nT,nParamTrys,pLo,pHi,uselogsampling);

% Evaluate FIM, calculate eigenvalues, and pick starting positions
% that appear to be far from plateaus based on the number of
% nonzero eigenvalues
%[Fs,eigvals] = deal(cell(nParamTrys,1));
nnzeigs = zeros(nParamTrys,1);
update = m.Update;
% Include constraint objs in the FIM calculation, if they are provided
if isfield(fitopts, 'ConstraintObj')
    obj_temp = vertcat(obj, fitopts.ConstraintObj{:});
else
    obj_temp = obj;
end
fitopts.Verbose = 0;
fprintf('Calculating FIMS for various starting parameter sets...\n')
pb = parforprogressbar(nParamTrys);
if nargout > 1
    Fs = cell(nParamTrys,1);
end
nargsout = nargout;
parfor i = 1:nParamTrys
    mi = update(p(:,i));
    if nargsout > 1
        [nnzeigs(i),Fs{i}] = calculateMetric(mi, con, obj_temp, fitopts, initopts);
    else
        nnzeigs(i) = calculateMetric(mi, con, obj_temp, fitopts, initopts);
    end
    pb.printbar(i)
end
[nnzeigs_sort,isort] = sort(nnzeigs,1,'descend');
fprintf('Best FIM has a rank score of %g\n', nnzeigs_sort(1))
k = p(:,isort(1:nParamSets));

end


function [val,F] = calculateMetric(m, con, obj, opts, initopts)

if initopts.UseHessian
    F = ObjectiveHessian(m, con, obj, opts);
else
    F = ObjectiveInformation(m, con, obj, opts);
end

threshold = initopts.Threshold;

switch initopts.Metric
    case 'SumOfConditionNumbers'
        eigvals = abs(eig(F));
        condnos = eigvals./max(eigvals);
        val = sum(condnos);
    case 'EigenvalueThresholdCount'
        eigvals = abs(eig(F));
        val = sum(eigvals > threshold);
    case 'ConditionNumberThresholdCount'
        eigvals = abs(eig(F));
        condnos = eigvals./max(eigvals);
        val = sum(condnos > threshold);
    case 'ObjectiveScaledEigenvalueThresholdCount'
        sim = SimulateSystem(m, con, obj, opts);
        int = [sim.int];
        for i = numel(int):-1:1
            err(i) = {obj(i).err(int(i))};
        end
        G = sum(cellfun(@(err)sum(err.^2), err));
        %G = ObjectiveValue(m, con, obj, opts);
        eigvals = abs(eig(F));
        val = sum(eigvals./G > threshold);
    case 'SumOfObjectiveScaledEigenvalues'
        G = ObjectiveValue(m, con, obj, opts);
        val = trace(abs(F))./G;
    case 'Random'
        val = rand;
    otherwise
        error('Unrecognized metric.')
end

end

