function monoval = ComputeMultiDegree(xval,degs)
% Make sure their size match
varnum = max(size(xval,1),size(xval,2));
monoval = 1;
for ind = 1:varnum
    monoval = monoval*xval(ind)^degs(ind);
end
end