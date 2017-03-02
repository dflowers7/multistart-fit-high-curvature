function p = randomlySampleParameters(nParams,nParamTrys,pLo,pHi)

% Fix 0 or Inf bounds
if any(pLo == 0)
    warning('Zero-valued lower bounds detected. Assuming a 1e-6 lower bound for log-space sampling.')
    pLo(pLo == 0) = 1e-6;
end
if any(isinf(pHi))
    warning('Inf-valued upper bounds detected. Assuming a 1e6 upper bound for log-space sampling.')
    pHi(isinf(pHi)) = 1e6;
end

deltap = bsxfun(@times, log10(pHi)-log10(pLo), rand(nParams, nParamTrys));

p = bsxfun(@plus, log10(pLo), deltap);
p = 10.^p;

end