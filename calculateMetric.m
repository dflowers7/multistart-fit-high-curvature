function rank = calculateMetric(F, G, metric, threshold)
% Input arguments
%   metric
%       'SumOfConditionNumbers'
%       'EigenvalueThresholdCount'
%       'ConditionNumberThresholdCount'
%       'ObjectiveScaledEigenvalueThresholdCount'
%       'SumOfObjectiveScaledEigenvalues'
%       'Objective
%       'Random'
if nargin == 0
    rank = {
        'SumOfConditionNumbers'
        'EigenvalueThresholdCount'
        'ConditionNumberThresholdCount'
        'ObjectiveScaledEigenvalueThresholdCount'
        'SumOfObjectiveScaledEigenvalues'
        'Random'
        };
    return
end

switch metric
    case 'SumOfConditionNumbers'
        eigvals = abs(eig(F));
        condnos = eigvals./max(eigvals);
        rank = sum(condnos);
    case 'EigenvalueThresholdCount'
        eigvals = abs(eig(F));
        rank = sum(eigvals > threshold);
    case 'ConditionNumberThresholdCount'
        eigvals = abs(eig(F));
        condnos = eigvals./max(eigvals);
        rank = sum(condnos > threshold);
    case 'ObjectiveScaledEigenvalueThresholdCount'
        %G = ObjectiveValue(m, con, obj, opts);
        eigvals = abs(eig(F));
        rank = sum(eigvals./G > threshold);
    case 'SumOfObjectiveScaledEigenvalues'
        rank = trace(abs(F))./G;
    case 'ObjectiveValue'
        rank = G;
    case 'Random'
        rank = rand;
    otherwise
        error('Unrecognized metric.')
end

end