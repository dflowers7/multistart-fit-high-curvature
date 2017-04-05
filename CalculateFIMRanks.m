function tab = CalculateFIMRanks(tab, rankmetrics, rankthresholds)

tab_old = tab;
nMetrics = numel(rankmetrics);
tab = cell(numel(tab_old),1);

Fs = {tab_old.Fs};
Gs = tab_old.Gs;

for i = 1:numel(tab_old)
    tab_i = tab_old(i);
    Fs = tab_i.Fs;
    Gs = tab_i.Gs;
    nFIMs_i = size(Fs,3);
    tab_i.rankstruct = struct('metric', cell(nMetrics,nFIMs_i), 'value', cell(nMetrics,nFIMs_i), 'threshold', cell(nMetrics,nFIMs_i));
    for j = nMetrics:-1:1
        rankvals_ij = cell(nFIMs_i,1);
        rankmetric_j = rankmetrics{j};
        rankthreshold_j = rankthresholds{j};
        for k = nFIMs_i:-1:1
            rankvals_ij{k} = calculateMetric(Fs(:,:,k), Gs(k), rankmetric_j, rankthreshold_j);
        end
        [tab_i.rankstruct(j,:).metric] = deal(rankmetric_j);
        [tab_i.rankstruct(j,:).value] = deal(rankvals_ij{:});
        [tab_i.rankstruct(j,:).threshold] = deal(rankthreshold_j);
    end
    tab{i} = tab_i;
end
tab = vertcat(tab{:});

end