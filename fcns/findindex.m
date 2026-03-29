function pos = findindex(a_f,v)
pos = 0;
for ind = 1:size(a_f,1)
    if isequal(v, a_f(ind, :))
        pos = ind;
        return;
    end
end
end