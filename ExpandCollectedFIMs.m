function tab = ExpandCollectedFIMs(tab)

metricNames = calculateMetric();
outputArgNames_all = {'ranks','ms','cons','Fs','rankstruct'};
outputArgNames = intersect(fieldnames(tab), outputArgNames_all);
inputArgNames = setdiff(fieldnames(tab), outputArgNames);

% Split jobs up into separate FIMs
njobs = size(tab,1);
tab_temp = cell(njobs,1);
for ji = 1:njobs
    tab_ji_unexpanded = tab(ji);
    nFIMs = tab_ji_unexpanded.nFIMs;
    for fi = nFIMs:-1:1
        % Input arguments are the same for each FIM
        for ai = 1:numel(inputArgNames)
            if ~isempty(inputArgNames{ai})
                tab_fi.(inputArgNames{ai}) = tab_ji_unexpanded.(inputArgNames{ai});
            end
        end
        % Output arguments vary for each FIM and are stored originally
        % as arrays. Split the arrays up into separate table entries.
        for ai = 1:numel(outputArgNames)
            if ~isempty(outputArgNames{ai})
                switch outputArgNames{ai}
                    % Experiments are indexed slightly differently than
                    % models and ranks
                    case 'cons'
                        tab_fi.(outputArgNames{ai}) = tab_ji_unexpanded.(outputArgNames{ai})(:,fi);
                    otherwise
                        tab_fi.(outputArgNames{ai}) = tab_ji_unexpanded.(outputArgNames{ai})(fi);
                end
            end
        end
        tab_ji(fi,1) = tab_fi;
    end
    tab_temp{ji} = tab_ji;
end
tab = vertcat(tab_temp{:});


end