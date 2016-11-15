function k = chooseInitialParametersByCurvature(m, con, obj, fitopts, initopts)
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

nParamSets = 1; % This function is currently used to select one parameter set only, but I'm leaving this here in case I want to add functionality

nParams = m.nk;
pLo = fitopts.LowerBound;
pHi = fitopts.UpperBound;

% Fix 0 or Inf bounds
if any(pLo == 0)
    warning('Zero-valued lower bounds detected. Assuming a 1e-6 lower bound for log-space sampling.')
    pLo(pLo == 0) = 1e-6;
end
if any(isinf(pHi))
    warning('Inf-valued upper bounds detected. Assuming a 1e6 upper bound for log-space sampling.')
    pHi(isinf(pHi)) = 1e6;
end

nParamTrys = 100*nParamSets;
p = randomlySampleParameters(nParams,nParamTrys,pLo,pHi);

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
parfor i = 1:nParamTrys
    mi = update(p(:,i));
    nnzeigs(i) = calculateMetric(mi, con, obj_temp, fitopts, initopts);
    pb.printbar(i)
end
[nnzeigs_sort,isort] = sort(nnzeigs,1,'descend');
fprintf('Best FIM has a rank score of %g\n', nnzeigs_sort(1))
k = p(:,isort(1:nParamSets));

end

function p = randomlySampleParameters(nParams,nParamTrys,pLo,pHi)

p = 10.^(bsxfun(@times, log10(pLo)+(log10(pHi)-log10(pLo)), rand(nParams, nParamTrys)));

end


function val = calculateMetric(m, con, obj, opts, initopts)

if initopts.UseHessian
    F = ObjectiveHessian(m, con, obj, opts);
else
    F = ObjectiveInformation(m, con, obj, opts);
end

threshold = initopts.Threshold;

switch initopts.Metric
    case 'SumOfConditionNumbers'
        eigvals = abs(eig(F));
        condnos = eigvals./eigvals(1);
        val = sum(condnos);
    case 'EigenvalueThresholdCount'
        eigvals = abs(eig(F));
        val = sum(eigvals > threshold);
    case 'ConditionNumberThresholdCount'
        eigvals = abs(eig(F));
        condnos = eigvals./eigvals(1);
        val = sum(condnos > threshold);
    case 'ObjectiveScaledEigenvalueThresholdCount'
        G = ObjectiveValue(m, con, obj, opts);
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

