function tab = CollectSampledFIMs(direcs, modelgenfun_arginnames, keepms, keepcons, inputArgFilter, inputargtables, JobFileOutputs, FIMskept, Gskept)
% Input arguments
%   inputArgFilter
%       .name   [ string ]
%       .values [ cell vector of values to load ]
%   inputargtables
%       cell array of input argument tables for each directory in direcs
%   JobFileOutputs
%       cell array of cell arrays of first six outputs of JobFiles()
%       function for each directory in direcs

if nargin < 9
    Gskept = [];
    if nargin < 8
        FIMskept = [];
        if nargin < 7
            JobFileOutputs = [];
            if nargin < 6
                inputargtables = [];
                if nargin < 5
                    inputArgFilter = [];
                    if nargin < 4
                        keepcons = [];
                        if nargin < 3
                            keepms = [];
                        end
                    end
                end
            end
        end
    end
end

if isempty(Gskept)
    Gskept = true;
end
if isempty(FIMskept)
    FIMskept = false;
end
if isempty(keepms)
    keepms = true;
end
if isempty(keepcons)
    keepcons = true;
end
direcs = cellstr(direcs);
if isempty(inputargtables)
    inputargtables = cell(numel(direcs),1);
end
if isempty(JobFileOutputs)
    JobFileOutputs = cell(numel(direcs),1);
end

% SampleFIMs_job(modelgenfun, mfile, rankfun, nFIMs, uselogsampling, useObj, useConstraints, varargin)
% Load inputs and outputs for each job
other_arginnames = {'nworkersperjob','fun','modelgenfun', 'mfile', 'rankfun', 'nFIMs', 'uselogsampling', 'useObj', ...
    'useConstraints', 'rngseed'};
inputArgNames = [other_arginnames(:)' modelgenfun_arginnames];
outputArgNames = {'ranks','ms','cons','Gs'};
if FIMskept
    outputArgNames = [outputArgNames {'Fs'}];
end
skippedInputArgs = {'nworkersperjob','fun'};
skippedOutputArgs = {};
if ~keepms
    skippedOutputArgs = [skippedOutputArgs {'ms'}];
end
if ~keepcons
    skippedOutputArgs = [skippedOutputArgs {'cons'}];
end

tab = cell(numel(direcs),1);
for di = 1:numel(direcs)
    % Load results for each job
    tab{di} = InputOutputTable(direcs{di}, inputArgNames, outputArgNames, skippedInputArgs, skippedOutputArgs, inputArgFilter, inputargtables{di}, JobFileOutputs{di});
end

tab = vertcat(tab{:});

end