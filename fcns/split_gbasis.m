function subbasis = split_gbasis(gb)

subbasis = [];
nGb = gb;
while size(nGb,2)>0
    tGb = nGb(1);
    index = 1;
    nGb(1) = [];
    subvar = symvar(tGb);
    resvar = symvar(nGb);
    comvar = intersect(subvar,resvar);
    while size(comvar,2)>0
        index = 1;
        while index <= size(nGb,2)
            nvar = symvar(nGb(index));
            intvar = intersect(subvar,nvar);
            if size(intvar,2)>0
                tGb = [tGb,nGb(index)];
                nGb(index) = [];
                subvar = union(subvar,nvar);
            else
                index = index + 1;
            end
        end
        resvar = symvar(nGb);
        comvar = intersect(subvar,resvar);
    end
    ncell = struct('Equations',tGb,'variables',subvar);
    subbasis = [subbasis,ncell];
end
end