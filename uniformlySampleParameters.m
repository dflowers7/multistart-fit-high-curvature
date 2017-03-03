function p = uniformlySampleParameters(nParams,nParamTrys,pLo,pHi,uselogsampling)

% Take log of parameters
pLo(uselogsampling) = log(pLo(uselogsampling));
pHi(uselogsampling) = log(pHi(uselogsampling));

% Fix Inf or -Inf bounds
if any(isinf(pHi))
    warning('Inf-valued upper bounds detected. Assuming a 1e6 upper bound.')
    pHi(isinf(pHi) & uselogsampling & sign(pHi)==1) = log(1e6);
    pHi(isinf(pHi) & ~uselogsampling & sign(pHi)==1) = 1e6;
    % Handles -Inf as a lower bound. You've probably got problems if this code is needed.
    pHi(isinf(pHi) & uselogsampling & sign(pHi)==-1) = log(1e-6);
    pHi(isinf(pHi) & ~uselogsampling & sign(pHi)==-1) = -1e6;
end
if any(isinf(pLo))
    warning('-Inf-valued lower bounds detected. Assuming a -1e6 lower bound for linear-basis parameters and a 1e-6 lower bound for log-basis parameters.')
    pLo(isinf(pLo) & uselogsampling & sign(pLo)==-1) = log(1e-6);
    pLo(isinf(pLo) & ~uselogsampling & sign(pLo)==-1) = -1e6;
    % Handles +Inf as a lower bound. You've probably got problems if this code is needed.
    pLo(isinf(pLo) & uselogsampling & sign(pLo)==1) = log(1e6);
    pLo(isinf(pLo) & ~uselogsampling & sign(pLo)==1) = 1e6; 
end

deltap = bsxfun(@times, pHi-pLo, rand(nParams, nParamTrys));

p = bsxfun(@plus, pLo, deltap);
p(uselogsampling,:) = exp(p(uselogsampling,:));

end