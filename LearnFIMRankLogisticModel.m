function rankmodel = LearnFIMRankLogisticModel(tab, numberParamSetsToKeep)
% Input arguments:
%   numberParamSetsToKeep
%       If set, will only use the top so many parameter sets to train the
%       model.

if isfield(tab, 'ms') && ~isfield(tab, 'ks')
    temp = vertcat(tab.ms);
    ks = [temp.k];
else
    ks = [tab.ks];
end

if isfield(tab, 'cons') && ~isfield(tab, 'ss')
    temp = [tab.cons];
    ss = {temp.ss};
    ncon = size(tab(1).cons, 1);
    ss = reshape(ss, ncon, numel(ss)/ncon);
    ss = cell2mat(ss);
else
    ss = [];
end

if isfield(tab, 'cons') && ~isfield(tab, 'qs')
    temp = [tab.cons];
    qs = {temp.qs};
    ncon = size(tab(1).cons,1);
    qs = reshape(qs, ncon, numel(qs)/ncon);
    qs = cell2mat(qs);
else
    qs = [];
end

if isfield(tab, 'cons') && ~isfield(tab, 'hs')
    temp = [tab.cons];
    hs = {temp.hs};
    ncon = size(tab(1).cons,1);
    hs = reshape(hs, ncon, numel(hs)/ncon);
    hs = cell2mat(hs);
else
    hs = [];
end

ps = {ks;ss;qs;hs};
ps = vertcat(ps{:}).';
logps = log(ps);

ranks = vertcat(tab.ranks);
[sortranks,sortranks_i] = sort(ranks,1,'descend');
i_keep = sortranks_i(1:numberParamSetsToKeep);
logps_keep = logps(i_keep,:);
ranks_keep = ranks(i_keep);

[B,dev,stats] = mnrfit(logps_keep, ranks_keep+1, 'Model', 'ordinal', 'Interactions', 'on');

rankmodel.predict = @(X)mnrval(B,log(X.'),stats,'Model','ordinal','Interactions','on');
rankmodel.B = B;
rankmodel.stats = stats;

end