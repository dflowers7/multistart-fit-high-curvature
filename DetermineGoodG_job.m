function [G,k,s,q,h] = DetermineGoodG_job(modelgenfun, varargin)
% [G,k,s,q,h] = DetermineGoodG_job(modelgenfun, rngseed, GRelTol, fitfun,
% useObj, useConstraints, initopts)

% Process inputs
p = inputParser();
p.addRequired('modelgenfun');
p.addOptional('rngseed', []);
p.addOptional('GRelTol', 0.05);
p.addOptional('fitfun', @FitObjective);
p.addOptional('useObj', []);
p.addOptional('useConstraints',[]);
p.addOptional('initopts',[]);
p.parse(modelgenfun, varargin{:});
res = p.Results;
modelgenfun = res.modelgenfun;
rngseed = res.rngseed;
GRelTol = res.GRelTol;
fitfun = res.fitfun;
useObj = res.useObj;
useConstraints = res.useConstraints;
initopts = res.initopts;

% Process empty inputs
if isempty(GRelTol)
    GRelTol = 0.01;
end

% Set random number seed, if provided
if ~isempty(rngseed)
    rngOrig = rng(rngseed);
end

% Generate model, experiment, objective, and options using the user-written
% routine
[m,con,obj,opts] = modelgenfun();

% Add on output function that stops fit once stagnation occurs
stagnationfun = StagnationOutputFunction(GRelTol);
if isfield(opts, 'OutputFcn')
    oldOutputFun = opts.OutputFcn;
else
    oldOutputFun = [];
end
if isempty(oldOutputFun)
    opts.OutputFcn = stagnationfun;
else
    opts.OutputFcn = @appendedStagnationOutputFunction;
end

    function stop = appendedStagnationOutputFunction(x, optimValues, state)
        stop = oldOutputFun(x, optimValues, state);
        if ~stop
            stop = stagnationfun(x, optimValues, state);
        end
    end

% Fit until stagnation
[m,con,G] = FitObjectiveFromHighCurvature(m, con, obj, opts, fitfun, useObj, useConstraints, initopts);

% Pull out parameter values
k = m.k;
s = [con.s];
q = {con.q};
h = {con.h};

% Reset rng to original state
if ~isempty(rngseed)
    rng(rngOrig);
end

end