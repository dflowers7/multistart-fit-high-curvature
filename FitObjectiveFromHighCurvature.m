function varargout = FitObjectiveFromHighCurvature(m, con, obj, opts, fitfun, useObj, useConstraints, initopts)

if nargin < 8
    initopts = [];
    if nargin < 7
        useConstraints = [];
        if nargin < 6
            useObj = [];
            if nargin < 5
                fitfun = [];
            end
        end
    end
end

if isempty(fitfun)
    fitfun = @FitObjective;
end
if isempty(useObj)
    useObj = true(size(obj));
end
if isempty(initopts)
    initopts.Metric = 'EigenvalueThresholdCount';
    initopts.UseHessian = false;
    initopts.Threshold = 1e-3;
end

hasConstraint = isfield(opts, 'ConstraintObj');

if hasConstraint && isempty(useConstraints)
    useConstraints = cell(size(opts.ConstraintObj));
    useConstraints(:) = {true(size(opts.ConstraintObj))};
end

obj_0 = obj;
opts_0 = opts;
objzero = objectiveZero();
obj_0(~useObj) = objzero;
if hasConstraint
    for i = 1:numel(opts.ConstraintObj)
        opts_0.ConstraintObj{i}(~useConstraints{i}) = objzero;
    end
end

% Randomly select 100 parameter sets, then pick the one with a FIM with the most
% nonzero eigenvalues
k0 = chooseInitialParametersByCurvature(m, con, obj_0, opts_0, initopts);

m = m.Update(k0);

varargout = cell(nargout,1);
[varargout{:}] = fitfun(m, con, obj, opts);

end