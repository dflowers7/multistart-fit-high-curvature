function submitOpts = defaultSubmitOpts(submitOpts, funname)

if nargin < 2
    funname = 'multistart_fit_high_curvature';
    if nargin < 1
        submitOpts = [];
    end
end

% Set defaults for cluster submission options
if isempty(submitOpts)
    submitOpts = struct;
end
default.name = funname;
default.proc = 4;
default.time = 24;
default.jvm = true;
submitOpts = mergestruct(default, submitOpts);

end