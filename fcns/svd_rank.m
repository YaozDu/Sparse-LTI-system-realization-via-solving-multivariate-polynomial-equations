function [rank,svars] = svd_rank(M,threshold,dec)
svars = svd(M);
rank = 0;
sp = svars(1);
for ind = 1:size(svars,1)
    sc = svars(ind);
    cdec = sc/sp;
    if sc>= threshold && cdec>=dec
        rank = rank+1;
        sp = sc;
    else
        break;
    end
end
end