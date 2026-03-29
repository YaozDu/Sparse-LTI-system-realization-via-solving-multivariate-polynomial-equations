function B = greedySieve(Mt,s,a_f,subsize,threshold,decay)
submtxsize = subsize(s); % This is the size of Ms(y), when s=1, we get size of M0(y)
varnum = size(a_f,2);
% L0 = a_f(1:submtxsize,:); % The candidates for B
str = struct('index','degree');

Msub = Mt(1:submtxsize,1:submtxsize); % The submatrix

% We can always add the first column, i.e. 1 into B. So we complete this
% first
Mp = Msub(:,1);
scanned = ones(submtxsize,1);
scanned(1) = 0;
str.index = 1;
str.degree = zeros(1,varnum);
B = str;

% Check remaining columns
for ind0 = 2:submtxsize
    if scanned(ind0) == 0
        continue % Skip terms that are eliminated
    end
    Col = Msub(:,ind0);
    Mc = [Mp,Col];
    Mc_col_size = size(Mc,2);
    [Mcrank,svar] = svd_rank(Mc,threshold,decay);
    if Mcrank == Mc_col_size
        Mp = Mc;
        scanned(ind0) = 0; % Add this term to B
        str.index = ind0;
        str.degree = a_f(ind0,:);
        B = [B,str];
    else
        % Time for elimination.
        degref = a_f(ind0,:);
        for ind2 = ind0:submtxsize % We remove unneeded terms
            cmp = (a_f(ind2,:)>=degref);
            if all(cmp)
                scanned(ind2) = 0;
            end
        end
    end
end
end