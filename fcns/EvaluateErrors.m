function cerr = EvaluateErrors(hs,zeval,a_f)
% Evaluate current value
Multidegrees = [];
for ind = 1:size(a_f,1)
    cdeg = a_f(ind,:);
    monomialvar = ComputeMultiDegree(zeval,cdeg);
    Multidegrees = [Multidegrees,monomialvar];
end
nvec = zeros(max(size(hs,1),size(hs,2)),1);
for ind = 1:max(size(hs,1),size(hs,2))
    hc = hs{ind};
    hcl = size(hc,1);
    nvec(ind) = Multidegrees(1:hcl)*hc;
end
cerr = (nvec'*nvec)^0.5;
end