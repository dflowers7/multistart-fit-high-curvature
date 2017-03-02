function [G,k] = DetermineGoodG_analyze(direc)

% Load outputs
outfiles = dir(fullfile(direc, 'output*.mat'));
outfiles = fullfile(direc, {outfiles.name});
for i = numel(outfiles):-1:1
    temp = load(outfiles{i});
    Gs(i) = temp.out{1};
    ks(:,i) = temp.out{2};
end

% Find the best objective value
[G,i_min] = min(Gs);
k = ks(:,i_min);

end