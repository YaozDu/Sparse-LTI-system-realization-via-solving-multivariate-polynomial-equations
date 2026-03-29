
function orders = getTermDegree(Term,x,varnum)
orders = zeros(1,varnum);
for ind = 1:varnum
    orders(ind) = polynomialDegree(Term,x(ind));
end
end