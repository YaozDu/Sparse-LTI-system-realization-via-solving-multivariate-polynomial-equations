% Given controller (Ap,Bp,Cp,0), and support S = (Asupp,Bsupp,Csupp,0),
% find realization (A,B,C,0) that satiesfies support S
clc,clear,close all

F = [0.2,0.6,0.6,0.2;0.6,0.8,0.8,0.6;0,0.4,0.6,0.2;0.6,0.4,0.4,0];
G = [0.8;0.2;0.4;0.6];
H = [0.4,0.8,0.4,0.6];

Ap = [-17.16,0,0,0;0,7.64,0,0;0,0,-0.48,0;0,0,0,-0.08];
Bp = [-9.2;5;-1.04;-0.76];
Cp = [13.68,-8.48,-0.6,0.28];

Asupp = [1,0,0,1;1,1,1,0;1,0,1,0;0,0,0,1];
Bsupp = [1;0;0;1];
Csupp = [1,1,0,0];

n = size(Asupp,1);

syms x [1,n*(n-1)] 
syms y [1,n*(n-1)] 
syms ytemp % This is for swapping variables

p = 2;
q = 1;
Ac = Ap;
Bc = Bp;
Cc = Cp;
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

Gb = gbasis(Eqns,x,'MonomialOrder','lexicographic');

Gblen = size(Gb,2);

%% Now, we want numbers
xsol = zeros(1,n*(n-1));
for ind = 1:n*(n-1)
    ind_elim = n*(n-1)-ind+1; 
    xelim = x(ind_elim);

    % Find intersection:
    Gblen = size(Gb,2);
    randgen = 1;
    cf = [];
    for ind2 = 1:Gblen
        cvars = symvar(Gb(ind2));
        if isequal(cvars,xelim)
            % Process the intersection variables
            fprintf('Eliminationg x(%d)\n',ind_elim);
            randgen = 0;
            cf = Gb(ind2);
            break
        end
    end

    % Eliminationg variables. randgen = 1, then we pick a variable. If
    % randgen = 0, we have to solve the variable
    if randgen
        % In this case, we can try assiging a value to the variable.
        cxval = 0;
        nGb = subs(Gb,xelim,cxval);
        nGb = gbasis(nGb,x,"MonomialOrder",'lexicographic');
        while isequal(nGb,sym(1))
            cxval = round(100*rand(1)-50);
            nGb = subs(Gb,xelim,cxval);
            nGb = gbasis(nGb,x,"MonomialOrder",'lexicographic');
        end
        % solution valid
        xsol(ind_elim) = cxval; % apply solution, eliminate the next one
        Gb = nGb;
    else
        % Now we need to solve the univariate polynomial
        cffs  = coeffs(cf);
        cffs = fliplr(cffs);
        cpMtx = (compan(cffs));
        if size(cpMtx,1) == 0 % In this case we have f = x(i)
            nGb = subs(Gb,xelim,0);
            Gb = gbasis(nGb,x,"MonomialOrder",'lexicographic'); % Update Gbasis
            xsol(ind_elim) = 0;
        else % We pick a valid solution
            xcand = eig(cpMtx);
            for ind3 = 1:size(xcand,1)
                cxval = xcand(ind3); % Try plugging this in
                nGb = subs(Gb,xelim,cxval);
                ceq = nGb(ind2);
                if size(symvar(ceq),2) == 0
                    ceq = double(ceq);
                    if abs(ceq)<=1e-8
                        nGb(ind2) = [];
                    end
                end
                if size(nGb,2)>0
                nGb = gbasis(nGb,x,"MonomialOrder",'lexicographic');
                if ~isequal(nGb,sym(1)) % solution valid
                    xsol(ind_elim) = cxval; % apply solution, eliminate the next one
                    Gb = nGb;
                    break;
                end
                else
                    xsol(ind_elim) = cxval;
                    break;
                end
            end
        end
    end
end

% Now, recover the solution:

T = eye(n);
iT = eye(n);

p = 2;
q = 1;
Ac = Ap;
Bc = Bp;
Cc = Cp;
for ind = 1:n*(n-1)
    [p,q];
    Ei = sym(eye(n));
    iEi = sym(eye(n));
    Ei(p,q) = xsol(ind);
    iEi(p,q) = -xsol(ind);
    T = T*Ei;
    iT = iEi*iT;
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

T = double(T);
iT = double(iT);

%% Display

Bc = double(Bc);
Cc = double(Cc);
Ac = double(Ac);

disp("Asupp=");
disp(Asupp);
disp('A=');
disp(Ac);

disp('Bsupp=');
disp(Bsupp);
disp('B=');
disp(Bc);

disp('Csupp=')
disp(Csupp);
disp('C=');
disp(Cc);

xsol