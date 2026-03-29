% Given controller (Ap,Bp,Cp,0), and support S = (Asupp,Bsupp,Csupp,0),
% find realization (A,B,C,0) that satiesfies support S
clc,clear,close all

F = [0.2,0.6,0.6,0.2;0.6,0.8,0.8,0.6;0,0.4,0.6,0.2;0.6,0.4,0.4,0];
G = [0.8;0.2;0.4;0.6];
H = [0.4,0.8,0.4,0.6];

Ap = [-17.16,0,0,0;0,7.64,0,0;0,0,-0.48,0;0,0,0,-0.08];
Bp = [-9.2;5;-1.04;-0.76];
Cp = [13.68,-8.48,-0.6,0.28];


Asupp = [1,0,0,1;1,1,0,0;0,0,1,1;1,0,1,1];
Bsupp = [1;0;0;1];
Csupp = [1,0,1,0];

[n,m] = size(Bp);
[p,~] = size(Cp);

A0 = Ap;
B0 = Bp;
C0 = Cp;

orders = perms([1:n]);

PermutationNum = size(orders,1);
feasibility = 0;
for ind = 1:PermutationNum
    Ad0 = zeros(n);
    Bd0 = zeros(n,m);
    Cd0 = zeros(p,n);

    corder = orders(ind,:);

    for ind2 = 1:n
        cindex = corder(ind2);
        Ad0(ind2,ind2) = A0(cindex,cindex);
        Bd0(ind2,:) = B0(cindex,:);
        Cd0(:,ind2) = C0(:,cindex);
    end


    syms x [1,n*(n-1)]
    syms y [1,n*(n-1)]
    syms ytemp % This is for swapping variables

    p = 2;
    q = 1;
    Ac = Ad0;
    Bc = Bd0;
    Cc = Cd0;
    for ind = 1:n*(n-1)
        Ei = sym(eye(n));
        iEi = sym(eye(n));
        Ei(p,q) = x(ind);
        iEi(p,q) = -x(ind);
        An = iEi*Ac*Ei;
        Bn = iEi*Bc;
        Cn = Cc*Ei;
        Ac = An;
        Bc = Bn;
        Cc = Cn;
        % Update p,q
        if p == n
            p = 1;
            q = q+1;
        else
            p = p+1;
            if q == p
                p = p+1;
            end
        end
    end

    Eqns = sym([]);
    for ind1 = 1:n
        for ind2 = 1:n
            if Asupp(ind1,ind2) == 0
                newconstraint = Ac(ind1,ind2);
                Eqns = [Eqns,newconstraint];
            end
        end
    end

    for ind1 = 1:n
        for ind2 = 1:size(Bsupp,2)
            if Bsupp(ind1,ind2) == 0
                newconstraint = Bc(ind1, ind2);
                Eqns = [Eqns,newconstraint];
            end
        end
        for ind2 = 1:size(Csupp,1)
            if Csupp(ind2,ind1) == 0
                newconstraint = Cc(ind2,ind1);
                Eqns = [Eqns,newconstraint];
            end
        end
    end

    Eqns = collect(Eqns,x);

    gb = gbasis(Eqns,x,'MonomialOrder','lexicographic');
    if ~isequal(gb,sym(1))
        feasibility = 1;
        break
    end
end

if feasibility == 0
    disp("No feasible Implementation!")
else
    disp("Complex Implementation Exists")
end