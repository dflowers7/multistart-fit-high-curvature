function [ranks,ms,cons,Gs,Fs] = CalculateFIMs_job(modelgenfun, mfile, rankfun, nFIMs, samplerfun, useObj, useConstraints, rngseed, varargin)
% [ranks,ms,cons,Gs,Fs] = CalculateFIMs_job(modelgenfun, mfile, rankfun, nFIMs, samplerfun, useObj, useConstraints, rngseed, varargin)
% Input arguments:
%   samplerfun
%       @(nparameters, nsamples, lowerbounds, upperbounds)
rngstate = rng(rngseed);

% Generate the model, experiments, objective, and options
[m,con,obj,opts] = modelgenfun(mfile, varargin{:});

if ~isfield(opts, 'UseParams')
    opts.UseParams = [];
end
if ~isfield(opts, 'UseSeeds')
    opts.UseSeeds = [];
end
if ~isfield(opts, 'UseInputControls')
    opts.UseInputControls = [];
end
if ~isfield(opts, 'UseDoseControls')
    opts.UseDoseControls = [];
end

% Fix the various options and inputs
[UseParams,UseSeeds,UseInputControls,UseDoseControls] = ...
    fixUses(m, con, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

nT = countFitParameters(UseParams, UseSeeds, UseInputControls, UseDoseControls);
TLowerBounds = collectFitBounds(opts.LowerBound, UseParams, UseSeeds, UseInputControls, UseDoseControls);
TUpperBounds = collectFitBounds(opts.UpperBound, UseParams, UseSeeds, UseInputControls, UseDoseControls);

Ts = samplerfun(nT, nFIMs, TLowerBounds, TUpperBounds);
nFIMs = size(Ts,2); % Recalculate in case the provided nFIMs isn't actually a number

% Filter out unused objective structs..
objzero = objectiveZero();
% ...in the objective function
if isempty(useObj)
    useObj = true(size(obj));
end
obj(~useObj) = objzero;

% ...in the constraints
hasConstraint = isfield(opts, 'ConstraintObj');
if hasConstraint
    % Default to using all constraints' obj structs when calculating the FIM
    if isempty(useConstraints)
        useConstraints = cell(size(opts.ConstraintObj));
        useConstraints(:) = {true(size(opts.ConstraintObj))};
    end
    
    % Replace unused constraint objectives with objectiveZero
    for i = 1:numel(opts.ConstraintObj)
        opts.ConstraintObj{i}(~useConstraints{i}) = objzero;
    end
    
    % Concatenate obj with constraint objs
    obj = vertcat(obj, opts.ConstraintObj{:});
end

ncon = numel(con);
ms = rmfield(m, setdiff(fieldnames(m), {'Name','k'}));
ms = repmat(ms, nFIMs, 1);
cons = rmfield(con, setdiff(fieldnames(con), {'Name','s','q','h'}));
cons = repmat(cons, 1, nFIMs);
[ms,cons] = distributeFitParameters(ms,cons,Ts,UseParams,UseSeeds,UseInputControls,UseDoseControls);

ranks = nan(nFIMs,1);
ks = cell(nFIMs,1);
[ss,qs,hs] = deal(cell(numel(con),nFIMs));
Gs = nan(nFIMs,1);
keepFIMs = nargout > 4;
if keepFIMs
    Fs = nan(nT,nT,nFIMs);
end

% Calculate the FIMs and their ranks
pb = parforprogressbar(nFIMs);
opts_ = opts;
opts_.Verbose = 0;
ncon = numel(con);
mupdate = m.Update;
conupdate = {con.Update};
parfor i = 1:nFIMs
    m_i = mupdate(ms(i).k);
    con_i = cell(ncon,1);
    for j = 1:ncon
        con_i{j} = conupdate{j}(cons(j,i).s, cons(j,i).q, cons(j,i).h);
    end
    con_i = vertcat(con_i{:});
    try
        sim = SimulateSensitivity(m_i, con_i, obj, opts_);
        int = reshape([sim.int], size(obj));
        int_cell = mat2cell(int, size(int,1), ones(size(int,2),1));
        F_i = ObjectiveInformation(m_i, con_i, obj, opts_, int_cell);
        err = cell(numel(int),1);
        for j = 1:numel(int)
            err(j) = {obj(j).err(int(j))};
        end
        G_i = sum(cellfun(@(err)sum(err.^2), err));
    catch ME
        fprintf(ME.getReport)
        fprintf('FIM calculation failed, continuing to next parameter set.\n')
        continue
    end
    rank_i = rankfun(F_i,G_i);
    if ~isempty(rank_i)
        ranks(i) = rank_i;
    end
    ks{i} = m_i.k;
    ss(:,i) = {con_i.s}';
    qs(:,i) = {con_i.q}';
    hs(:,i) = {con_i.h}';
    Gs(i) = G_i;
    %[ss{:,i}, qs{:,i}, hs{:,i}] = deal(con_i.s, con_i.q, con_i.h);
    if keepFIMs
        Fs(:,:,i) = F_i;
    end
    pb.printbar(i)
end

% Set rng state back to its original value
rng(rngstate)

end