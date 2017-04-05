function tab = CompareFitsFromSampledFIMs(direc, modelgenfun_inputargnames)

modelgenfun_inputargnames = modelgenfun_inputargnames(:)';

[outfiles, outfiles_jis, infiles, infiles_jis, cmdfiles, cmdfiles_jis] = JobFiles(direc);

inputArgNames = [{'nworkersperjob','fun','modelgenfun', 'mfile', 'ks', 'ss', 'qs', 'hs', 'Gstop'}, modelgenfun_inputargnames];
outputArgNames = {'mbest', 'conbest', 'Gbest', 'Dbest'};
skippedInputArgs = {'nworkersperjob','fun','ks','ss','qs','hs'};
skippedOutputArgs = [];
tab = InputOutputTable(direc, inputArgNames, outputArgNames, skippedInputArgs, skippedOutputArgs);

end