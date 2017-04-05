function varargout = FitObjectiveFromHighCurvature(m, con, obj, opts, fitfun, useObj, useConstraints, initopts)
% varargout = FitObjectiveFromHighCurvature(m, con, obj, opts, fitfun, useObj, useConstraints, initopts)
% Returns same outputs as FitObjective, plus a fifth output returning the
% FIM of the starting parameter set of the optimization

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

% Fix FitObjective opts
[m,con,obj,opts] = FixFitObjectiveOpts(m,con,obj,opts);
opts.Verbose = opts.Verbose + 1; % Hack to keep Verbose from being decreased twice. Forces Verbose to be on, though.

hasConstraint = isfield(opts, 'ConstraintObj');

% Default to appending all constraints' functions when calculating the FIM
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
if nargout > 4
    [k0,s0,q0,h0,F] = chooseInitialParametersByCurvature(m, con, obj_0, opts_0, initopts);
else
    [k0,s0,q0,h0] = chooseInitialParametersByCurvature(m, con, obj_0, opts_0, initopts);
end

m = m.Update(k0);
for i = 1:numel(con)
    con(i) = con(i).Update(s0{i},q0{i},h0{i});
end

% Feig = eig(F);
% eigindex = sum(Feig > initopts.Threshold);
% %eigindex = 63;
% 
% % Set up constraint function that prevents optimization from reducing the
% % number of significant FIM eigenvalues
% [intfun,objfun] = GenerateFIMEigenvalueFunction(eigindex, true);
% opts.ConstraintObj = {obj};
% opts.ConstraintVal = -log(initopts.Threshold);
% opts.ConstraintIntegrateFunction = intfun;
% opts.ConstraintReductionFunction = objfun;


varargout = cell(nargout,1);
[varargout{1:min(nargout,4)}] = fitfun(m, con, obj, opts);
if nargout > 4
    varargout{5} = F;
end

end