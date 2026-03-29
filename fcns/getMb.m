function Mb = getMb(Mt,B)
Bsize = size(B,2);
Mb = zeros(Bsize);
for ind1 = 1:Bsize
    for ind2 = 1:Bsize
        xpos = B(ind1).index;
        ypos = B(ind2).index;
        Mb(ind1, ind2) = Mt(xpos, ypos);
    end
end
end