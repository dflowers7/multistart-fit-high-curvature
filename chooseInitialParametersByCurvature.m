function [k,s,q,h,Fs,p_k,p_s,p_q,p_h] = chooseInitialParametersByCurvature(m, con, obj, fitopts, initopts, uselogsampling)
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

nParamSetsToKeep = 1; % This function is currently used to select one parameter set only, but I'm leaving this here in case I want to add functionality

% Calculate the number of fit parameters
nTk = sum(fitopts.UseParams);
nTs_icon = sum(fitopts.UseSeeds, 1);
nTs = sum(nTs_icon);
nTq_icon = cellfun(@sum, fitopts.UseInputControls);
nTq = sum(nTq_icon);
nTh_icon = cellfun(@sum, fitopts.UseDoseControls);
nTh = sum(nTh_icon);

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
    nParamTrys = 100*nParamSetsToKeep;
end
T = uniformlySampleParameters(nT,nParamTrys,pLo,pHi,uselogsampling);
p_k = repmat(m.k, 1, nParamTrys);
p_k(fitopts.UseParams,:) = T(1:nTk,:);
[p_s, p_q, p_h] = deal(cell(numel(con),nParamTrys));
for i = 1:numel(con)
    p_s(i,:) = {con(i).s};
    p_q(i,:) = {con(i).q};
    p_h(i,:) = {con(i).h};
    
    UseSeeds_i = fitopts.UseSeeds(:,i);
    startind_s = nTk + sum(nTs_icon(1:i-1)) + 1;
    endind_s = startind_s + nTs_icon(i) - 1;
    
    UseInputControls_i = fitopts.UseInputControls{i};
    startind_q = nTk + nTs + sum(nTq_icon(1:i-1)) + 1;
    endind_q = startind_q + nTq_icon(i) - 1;
    
    UseDoseControls_i = fitopts.UseDoseControls{i};
    startind_h = nTk + nTs + nTq + sum(nTh_icon(1:i-1)) + 1;
    endind_h = startind_h + nTh_icon(i) - 1;
    for j = 1:nParamTrys
        % Assign seeds per experiment
        p_s{i,j}(UseSeeds_i) = T(startind_s:endind_s,j);
        
        % Assign input controls per experiment
        p_q{i,j}(UseInputControls_i) = T(startind_q:endind_q,j);
        
        % Assign dose controls per experiment
        p_h{i,j}(UseDoseControls_i) = T(startind_h:endind_h,j);
    end
end

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
parfor iTry = 1:nParamTrys
    update_i = update;
    m_iTry = update_i(p_k(:,iTry));
    con_iTry = con;
    for icon = 1:numel(con_iTry)
        con_iTry(icon) = con_iTry(icon).Update(p_s{icon,iTry}, p_q{icon,iTry}, p_h{icon,iTry});
    end
    if nargsout > 1
        [nnzeigs(iTry),Fs{iTry}] = calculateMetric(m_iTry, con_iTry, obj_temp, fitopts, initopts);
    else
        nnzeigs(iTry) = calculateMetric(m_iTry, con_iTry, obj_temp, fitopts, initopts);
    end
    pb.printbar(iTry)
end
[nnzeigs_sort,isort] = sort(nnzeigs,1,'descend');
fprintf('Best FIM has a rank score of %g\n', nnzeigs_sort(1))
k = p_k(:,isort(1:nParamSetsToKeep));
s = p_s(:,isort(1:nParamSetsToKeep));
q = p_q(:,isort(1:nParamSetsToKeep));
h = p_h(:,isort(1:nParamSetsToKeep));

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

