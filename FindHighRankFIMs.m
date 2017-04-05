function [ranks,ms,cons,Gs,Fs] = FindHighRankFIMs(cluster, nFIMsamples, napproxsamples, classifier, njobs, rankmetric, rankthreshold, modelgenfun, mfile, submitOpts, maxWait, uselogsampling, useObj, useConstraints, varargin)

input = [{cluster, nFIMsamples, napproxsamples, classifier, njobs, rankmetric, rankthreshold, modelgenfun, mfile, submitOpts, maxWait, uselogsampling, useObj, useConstraints} varargin(:)'];
save input1.mat input

classifier = validatestring(classifier, {'svm','naive-bayes'});

highestrank = 0;
nstagnatedruns = 0;
rankisimproving = true;

[m,con,obj,opts] = modelgenfun(mfile, varargin{:});

% Sample parameters uniformly for first round
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
[UseParams,UseSeeds,UseInputControls,UseDoseControls] = ...
    fixUses(m, con, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
nT = countFitParameters(UseParams, UseSeeds, UseInputControls, UseDoseControls);
TLowerBounds = collectFitBounds(opts.LowerBound, UseParams, UseSeeds, UseInputControls, UseDoseControls);
TUpperBounds = collectFitBounds(opts.UpperBound, UseParams, UseSeeds, UseInputControls, UseDoseControls);
Ts_nextround = uniformlySampleParameters(nT, nFIMsamples, TLowerBounds, TUpperBounds, uselogsampling);

Fs = [];
Gs = [];
Ts = zeros(nT,0);
Feigs = {};
ranks = zeros(0,1);
rank_thresh = 0;
last_rank_thresh = 0;
nstagnatedruns = 0;

while rankisimproving
    
    [Fs_temp,Gs_temp,~,~,Ts_nextround] = CalculateFIMs(cluster, modelgenfun, mfile, njobs, Ts_nextround, submitOpts, maxWait, useObj, useConstraints, varargin{:});
    Fs = cat(3, Fs, Fs_temp);
    Gs = vertcat(Gs, Gs_temp);
    Ts = [Ts Ts_nextround];
    
    Fs_temp_cell = mat2cell(Fs_temp,size(Fs_temp,1),size(Fs_temp,2),ones(size(Fs_temp,3),1));
    Fs_temp_cell = Fs_temp_cell(:);
    Gs_temp_cell = mat2cell(Gs_temp, ones(size(Gs_temp)));
    Feigs_temp = cellfun(@eig, Fs_temp_cell, 'UniformOutput', false);
    ranks_temp = cellfun(@(F,G)calculateMetric(F,G,rankmetric,rankthreshold), Fs_temp_cell, Gs_temp_cell);
    
    Feigs = [Feigs; Feigs_temp];
    ranks = [ranks; ranks_temp];
    
    ranks_sort = sort(ranks, 1, 'descend');
    
    % Get upper 95% value on sampled ranks, and use that as the threshold
    % for determining a "high" rank value
    %rank_thresh = ranks_sort(ceil(0.05*size(Ts,2)));
    % ...or select the nTth sorted rank, so that there are at least nT
    % samples with that rank or above
    last_rank_thresh = rank_thresh;
    rank_thresh = ranks_sort(nT);
    
    fraction_newhighranks = sum(ranks_temp >= rank_thresh)/nFIMsamples;
    
    maxrank = ranks_sort(1);
%     if maxrank > highestrank
%         highestrank = maxrank;
%         nstagnatedruns = 0;
%     else
%         nstagnatedruns = nstagnatedruns + 1;
%     end
    if rank_thresh > last_rank_thresh
        nstagnatedruns = 0;
    else
        nstagnatedruns = nstagnatedruns + 1;
    end
    rankisimproving = nstagnatedruns < 3;
    fprintf('Highest rank so far: %g\nMean rank: %g\nMean rank among new samples: %g\nHigh rank threshold: %g\nFraction of new samples that were high rank: %g\n', maxrank, mean(ranks_sort), mean(ranks_temp), rank_thresh, fraction_newhighranks)
    
    logTs = log10(Ts);
    
    switch classifier
        case 'svm'
            mdl = fitcsvm(logTs.', ranks>=rank_thresh, 'KernelFunction', 'gaussian', 'BoxConstraint', Inf, 'KernelScale', 'auto');
        case 'naive-bayes'    
%             rank_counts = accumarray(ranks+1,1);
%             rank_counts_to_eliminate = find(rank_counts==1)-1; % Cannot train a naive bayes classifier on a class with only one sample
%             elim = ismember(ranks, rank_counts_to_eliminate);
%             mdl = fitcnb(logTs(:,~elim).', ranks(~elim), 'DistributionNames', 'normal');

            mdl = fitcnb(logTs.', ranks>=rank_thresh, 'DistributionNames', 'normal');
        otherwise
            error('FindHighRankFIMs:UnrecognizedClassifier', 'Unrecognized classifier ''%s''.', classifier)
    end
    
    %mdl = fitcnb(logks(:,~elim).', ranks(~elim), 'DistributionNames', 'kernel'); % 'kernel' takes too long with 1000000 approximated samples
    
    %mdl = fitrgp(logks.', ranks, 'Sigma', 0.1, 'SigmaLowerBound', 0.00000001, 'KernelParameters', [mean(pdist(logks.'))*3;0.1]);
    
    %ishighrank = ranks >= rank_thresh;
    %ctree = fitctree(logks', ishighrank, 'PredictorNames', {m.Parameters.Name}');
    
    %mdl = fitctree(logks.', ranks);
    
    haveenoughsamples = false;
    Ts_nextround = nan(nT,nFIMsamples);
    nHighRankPointsToKeep = floor(nFIMsamples/2);
    nClosePointsToKeep = nFIMsamples-nHighRankPointsToKeep;
    Ts_nextround_highranks = nan(nT,nHighRankPointsToKeep);
    %Ts_nextround_almosthighranks = nan(nT,nAlmostHighRankPointsToKeep);
    Ts_nextround_close = nan(nT,nClosePointsToKeep);
    if isscalar(uselogsampling)
        uselogsampling_vec = repmat(uselogsampling, numel(TUpperBounds), 1);
    else
        uselogsampling_vec = uselogsampling;
    end
    temp_u = TUpperBounds;
    temp_u(uselogsampling_vec) = log10(temp_u(uselogsampling_vec));
    temp_l = TLowerBounds;
    temp_l(uselogsampling_vec) = log10(temp_l(uselogsampling_vec));
    logdT = (temp_u - temp_l)/2;
    lbounds_adjusted = TLowerBounds;
    ubounds_adjusted = TUpperBounds;
    while ~haveenoughsamples
        
        switch classifier
            case 'svm'
                
                % Sample parameters uniformly within the bounds
                Ts_test = uniformlySampleParameters(nT, napproxsamples, lbounds_adjusted, ubounds_adjusted, uselogsampling);
                logTs_test = Ts_test;
                logTs_test(uselogsampling_vec,:) = log10(logTs_test(uselogsampling_vec,:));
                
                % Use the prediction model to predict which parameter sets
                % will have high ranks
                [ranks_test, scores_test] = predict(mdl, logTs_test.');
                ishighrank_predicted_test = ranks_test > 0;
                
%                 [~, sortscores_test_i] = sort(scores_test(~ishighrank_predicted_test,2), 1, 'descend');
%                 nfound_almost = numel(sortscores_test_i);
%                 nleft_almost = sum(isnan(Ts_nextround_almosthighranks(1,:)));
%                 almosthighrank_i = find(~ishighrank_predicted_test);
%                 almosthighrank_i = almosthighrank_i(sortscores_test_i);
%                 almosthighrank_i_chosen = almosthighrank_i(sortscores_test_i(1:min(nleft_almost,nfound_almost)));
                
                % Keep parameter sets close to the SVM boundary, i.e.,
                % those with a score less than 0.5
                closescorethresh = 0.5;
                [sortabsscores, sortabsscores_test_i] = sort(abs(scores_test(:,2)));
                sortabsscores_test_i_chosen = sortabsscores_test_i(sortabsscores < closescorethresh);
                nfound_close = numel(sortabsscores_test_i_chosen);
                nleft_close = sum(isnan(Ts_nextround_close(1,:)));
                
                % Store parameter sets predicted to have a high rank
                i1_high = find(isnan(Ts_nextround_highranks(1,:)), 1);
                if ~isempty(i1_high)
                    % If is empty, then we've already filled up the matrix
                    % with entries
                    highrank_i = find(ishighrank_predicted_test, nHighRankPointsToKeep-i1_high+1);
                    i2_high = i1_high + numel(highrank_i) - 1;
                    Ts_nextround_highranks(:,i1_high:i2_high) = Ts_test(:,highrank_i);
                end
                
                % Store parameter sets predicted to be close to the SVM
                % boundary
                i1_close = find(isnan(Ts_nextround_close(1,:)), 1);
                if ~isempty(i1_close)
                    i2_close = i1_close + min(nfound_close,nleft_close) - 1;
                    Ts_nextround_close(:,i1_close:i2_close) = Ts_test(:,sortabsscores_test_i_chosen(1:min(nfound_close,nleft_close)));
                end
                
                % Once we have enough high-rank parameter sets and enough
                % close-to-boundary parameter sets, we're done with this
                % round
                haveenoughhighranksamples = i2_high >= nHighRankPointsToKeep;
                haveenoughclosesamples = i2_close >= nClosePointsToKeep;
                
                if ~haveenoughhighranksamples
%                     % If we don't have enough samples yet, adjust bounds to focus on points that are close to
%                     % the high-rank parameter set space
%                     [maxscore, maxscore_i] = max(scores_test(:,2));
%                     newlogmean = Ts_test(:,maxscore_i);
%                     newlogmean(uselogsampling_vec,:) = log10(newlogmean(uselogsampling_vec,:));
%                     logdT = logdT./2;
%                     lbounds_adjusted = newlogmean-logdT;
%                     lbounds_adjusted(uselogsampling_vec) = 10.^lbounds_adjusted(uselogsampling_vec);
%                     ubounds_adjusted = newlogmean+logdT;
%                     ubounds_adjusted(uselogsampling_vec) = 10.^ubounds_adjusted(uselogsampling_vec);
%                     lbounds_adjusted(lbounds_adjusted < TLowerBounds) = TLowerBounds(lbounds_adjusted < TLowerBounds);
%                     ubounds_adjusted(ubounds_adjusted > TUpperBounds) = TUpperBounds(ubounds_adjusted > TUpperBounds);

                    % If we don't have enough high-rank samples yet,
                    % randomly pair high-rank parameter sets with other
                    % parameter sets and sample for more points along the
                    % lines connecting them
                    linelogTs = [logTs logTs_test];
                    ishighrank_line = [ranks>=rank_thresh; ishighrank_predicted_test];
                    nhighrank = sum(ishighrank_line);
                    % If no high rank parameter sets found, just
                    % continue to the next round and keep sampling
                    % until we find at least one high-rank parameter set
                    if nhighrank > 0
                        while ~haveenoughhighranksamples
                            highrank_i = find(ishighrank_line);
                            nhighrank = numel(highrank_i);
                            nallrank = numel(ishighrank_predicted_test);
                            npairs = nhighrank;
                            highpairs = [highrank_i(randi(nhighrank, npairs, 1)) randi(nallrank, npairs, 1)];
                            nlinesamples = 10;
                            minfactor = 0.2; % Prevents selecting essentially the same point as the first point
                            maxfactor = 0.8; % Prevents selecting essentially the same point as the second point
                            linescalefactors = minfactor + (maxfactor-minfactor).*rand(npairs, nlinesamples);
                            linelogTs_new = nan(nT, npairs*nlinesamples);
                            for li = 1:npairs*nlinesamples
                                [pair_i, lscf_i] = ind2sub([npairs,nlinesamples], li);
                                logT1 = linelogTs(:,highpairs(pair_i,1));
                                logT2 = linelogTs(:,highpairs(pair_i,2));
                                linelogTs_new(:,li) = logT1 + (logT2-logT1)*linescalefactors(li);
                            end
                            % Shuffle the sampled parameters so that we don't
                            % tend to keep parameter sets drawn from the same
                            % line
                            linelogTs_new = linelogTs_new(:,randperm(npairs*nlinesamples));
                            ishighrank_line_new = predict(mdl, linelogTs_new.');
                            ishighrank_line_new = ishighrank_line_new > 0;
                            i1_high = i2_high+1;
                            i2_high = min(i1_high+sum(ishighrank_line_new)-1, nHighRankPointsToKeep);
                            nlinepointstokeep = i2_high-i1_high+1;
                            Ts_nextround_highranks(:,i1_high:i2_high) = 10.^linelogTs_new(:,find(ishighrank_line_new, nlinepointstokeep));
                            
                            haveenoughhighranksamples = i2_high >= nHighRankPointsToKeep;
                            
                            linelogTs = [linelogTs linelogTs_new];
                            ishighrank_line = [ishighrank_line; ishighrank_line_new];
                        end
                    end
                end
                
                if ~haveenoughclosesamples
                    % If not enough samples have been found that are close
                    % to the SVM border, randomly pair high-rank parameter
                    % sets with low-rank parameter sets and sample along
                    % the connecting lines
                    linelogTs = [logTs logTs_test];
                    ishighrank_line = [ranks>=rank_thresh; ishighrank_predicted_test];
                    nhighrank = sum(ishighrank_line);
                    nlowrank = sum(~ishighrank_line);
%                     close_i = sortabsscores_test_i_chosen;
%                     far_i = setdiff((1:size(linelogTs,2)).', close_i);
%                     nclose = numel(close_i);
%                     nfar = numel(far_i);
                    % If no parameter sets found close to the SVM border, just
                    % continue to the next round and keep sampling
                    % until we find at least one close parameter set
                    if nhighrank > 0
                        while ~haveenoughclosesamples
                            high_i = find(ishighrank_line);
                            low_i = find(~ishighrank_line);
                            npairs = nhighrank;
                            highlowrankpairs = [high_i(randi(nhighrank, npairs, 1)) low_i(randi(nlowrank, npairs, 1))];
                            nlinesamples = 10;
                            minfactor = 0.2;
                            maxfactor = 0.8;
                            linescalefactors = minfactor + (maxfactor-minfactor).*rand(npairs, nlinesamples);
                            linelogTs_new = nan(nT, npairs*nlinesamples);
                            for li = 1:npairs*nlinesamples
                                [pair_i, lscf_i] = ind2sub([npairs,nlinesamples], li);
                                logT1 = linelogTs(:,highlowrankpairs(pair_i,1));
                                logT2 = linelogTs(:,highlowrankpairs(pair_i,2));
                                linelogTs_new(:,li) = logT1 + (logT2-logT1)*linescalefactors(li);
                            end
                            % Shuffle the sampled parameters so that we don't
                            % tend to keep parameter sets drawn from the same
                            % line
                            linelogTs_new = linelogTs_new(:,randperm(npairs*nlinesamples));
                            [ishighrank_line_new,linescores_new] = predict(mdl, linelogTs_new.');
                            ishighrank_line_new = ishighrank_line_new > 0;
                            linescores_new = linescores_new(:,2);
                            close_i_new = find(abs(linescores_new) < closescorethresh);
                            i1_close = i2_close+1;
                            i2_close = min(i1_close+numel(close_i_new)-1, nClosePointsToKeep);
                            nlinepointstokeep = i2_close-i1_close+1;
                            Ts_nextround_close(:,i1_close:i2_close) = 10.^linelogTs_new(:,close_i_new(1:nlinepointstokeep));
                            
                            haveenoughclosesamples = i2_close >= nClosePointsToKeep;
                            
                            ishighrank_line = [ishighrank_line; ishighrank_line_new];
                            nhighrank = sum(ishighrank_line);
                            linelogTs = [linelogTs linelogTs_new];
                            
                        end
                    end
                end
                
                haveenoughsamples =  haveenoughhighranksamples && haveenoughclosesamples;
                
                if haveenoughsamples
                    Ts_nextround = [Ts_nextround_highranks Ts_nextround_close];
                end
            case 'naive-bayes'
                % Sample parameters from the distribution of the
                % high-rank (second) class so that more parameter sets are
                % positive hits
                mu_sd = mdl.DistributionParameters(end,:);
                mu_sd = [mu_sd{:}];
                mu = mu_sd(1,:);
                sd = mu_sd(2,:);
                zscores = randn(napproxsamples, nT);
                logTs_test = bsxfun(@plus, bsxfun(@times, zscores, sd), mu);
                Ts_test = exp(logTs_test).';
                ishighrank_predicted_test = predict(mdl, logTs_test);
                
                firstnan_i = find(isnan(Ts_nextround(1,:)), 1);
                highrank_i = find(ishighrank_predicted_test, nFIMsamples-firstnan_i+1);
                i1 = firstnan_i;
                i2 = firstnan_i + numel(highrank_i) - 1;
                Ts_nextround(:,i1:i2) = Ts_test(:,highrank_i);
                
                haveenoughsamples = i2 >= nFIMsamples;
        end

        
    end
    
end

[ms,cons] = distributeFitParameters([], [], Ts, UseParams, UseSeeds, UseInputControls, UseDoseControls);
[ms.Name] = deal(m.Name);
for i = 1:numel(con)
    [cons(i,:).Name] = deal(con(i).Name);
end