function [Chi,Soln_cand,commutability_mtx,sumerr,avg_err]=getChi(Px,Mb)
varnum = size(Px,3);
Bsize = size(Px,1);
Chi = zeros(Bsize,Bsize,varnum);
Soln_cand = zeros(Bsize,varnum);
commutability_mtx = zeros(varnum);
sumerr = 0;
maxComErr = 0;
avg_err = 0;

for ind = 1:varnum
    Chi(:,:,ind) = Mb^(-1)*Px(:,:,ind);
    Soln_cand(:,ind) = eig(Chi(:,:,ind));
end

if varnum>=2
    for ind1 = 1:varnum
        for ind2 = ind1+1:varnum
            eChi = Chi(:,:,ind1)*Chi(:,:,ind2)-Chi(:,:,ind2)*Chi(:,:,ind1);
            neChi = norm(eChi);
            if neChi>maxComErr
                maxComErr = neChi;
            end
            commutability_mtx(ind1,ind2) = neChi;
            sumerr = sumerr+neChi;
        end
    end
    avg_err = sumerr/(varnum*(varnum-1)/2);
end
end