function Px = getPx(Mt,B,a_f)
varnum = size(a_f,2);
Bsize = size(B,2);
Px = zeros(Bsize,Bsize,varnum);
for ind0 = 1:varnum
    dxi = zeros(1,varnum);
    dxi(ind0) = 1;
    for ind1 = 1:Bsize
        % xdeg = B(ind1).degree
        xpos = B(ind1).index;
        for ind2 = 1:Bsize
            % ydegbase = B(ind2).degree;
            ydeg = B(ind2).degree+dxi;
            ypos = findindex(a_f,ydeg);
            Px(ind1, ind2,ind0) = Mt(xpos, ypos);
        end
    end
end
end