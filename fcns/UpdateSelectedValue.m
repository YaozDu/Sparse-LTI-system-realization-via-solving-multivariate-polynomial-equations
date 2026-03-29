function nsel = UpdateSelectedValue(csel,varnum,Bsize)
for j = varnum:-1:1
    if csel(j)==Bsize && j>1
        csel(j) = 1;
    else
        csel(j) = csel(j)+1;
        break
    end
end
nsel = csel;
end