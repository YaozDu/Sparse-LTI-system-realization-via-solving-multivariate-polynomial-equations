% Given controller (Ap,Bp,Cp,0), and support S = (Asupp,Bsupp,Csupp,0),
% find realization (A,B,C,0) that satiesfies support S.
clc,clear,close all

F = [0.2,0.6,0.6,0.2;0.6,0.8,0.8,0.6;0,0.4,0.6,0.2;0.6,0.4,0.4,0];
G = [0.8;0.2;0.4;0.6];
H = [0.4,0.8,0.4,0.6];

Ap = [-17.16,0,0,0;0,7.64,0,0;0,0,-0.48,0;0,0,0,-0.08];
Bp = [-9.2;5;-1.04;-0.76];
Cp = [13.68,-8.48,-0.6,0.28];

Asupp = [1,1,0,1;1,1,0,0;0,0,1,1;0,0,1,1];
Bsupp = [1;0;0;1];
Csupp = [1,0,1,0];

setrng = 1121; % This ensures the right value is generated

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

gb = gbasis(Eqns,x,'MonomialOrder','lexicographic');
subbasis = split_gbasis(gb);

xsize = size(x,2);
x_sel = inf*zeros(1,xsize);

newsubbasis = [];

for index = 1:size(subbasis,2)
Gb0 = subbasis(index).Equations;
Gbvars = subbasis(index).variables;
Gbsize = size(Gb0,2);

[zerodim,xd] = detectIdealDim(Gb0,Gbvars);

Gbc = Gb0;
while ~zerodim
    % Substitude
    xsubs = xd(1);
    for le = 1:xsize
        if xsubs == x(le)
            break
        end
    end
    % Then we are assigning a number for x(le)

    % As always, try 0 first.
    xcand = 0;
    % xcand = round(20*rand(1)-10);
    nGb = subs(Gbc,xsubs,xcand);
    nGb = gbasis(nGb,Gbvars,'MonomialOrder','lexicographic');

    while isequal(nGb,sym(1))
        rng(setrng);
        xcand = round(400*rand(1)-200)/10;
        nGb = subs(Gbc,xsubs,xcand);
        nGb = gbasis(nGb,Gbvars,'MonomialOrder','lexicographic');
    end
    x_sel(le) = xcand;
    Gbc = nGb;
    Gbvars = setdiff(Gbvars,xsubs); % Kick out eliminated variable

    % Re-verify
    [zerodim,xd] = detectIdealDim(Gbc,Gbvars);
end

% Gbc is a zero dimensional ideal, we can try to split it.
newstr = split_gbasis(Gbc);
newsubbasis = [newsubbasis,newstr];
end

subbasis = newsubbasis;

% xsize = size(x,2);
% Gbc = subs(gb,x(1),-1.1);
% x_sel = inf*zeros(1,xsize);
% x_sel(1) = -1.1;
% Gbc = gbasis(Gbc,x,"MonomialOrder",'lexicographic');
% 
% % separate the grobner basis
% subbasis = split_gbasis(Gbc);

fragments = max(size(subbasis,2),size(subbasis,1));

threshold = 1e-7;
decay = 1e-3;

xsol = x_sel;
for iterations = 1:fragments
    Gbc = subbasis(iterations).Equations;
    redvars = subbasis(iterations).variables;
    varnum = size(redvars,2);
    backindex = zeros(1,varnum);
    % Find variables that z(j) correspond to.
    for i = 1:varnum
        var = redvars(i);
        for position = 1:size(x,2)
            if isequal(x(position), var)
                backindex(i) = position;
                break
            end
        end
    end

    syms z [1,varnum]
    Zeqs = subs(Gbc,redvars,z);

    degs = polynomialDegree(Zeqs);
    maxdeg= max(degs);
    d = ceil(maxdeg/2);
    dj = ceil(degs/2);

    td = 0;
    flag = 1;
    while flag == 1
        t = d+td; 
        [a_f,subsize] = multi_index_grlex(varnum,2*t);

        % Generate a base h(x):
        hs = cell(size(Zeqs,2),1);

        for ind = 1:size(Zeqs,2)
            h_i = Zeqs(ind);
            % The vector storing variables of h_i:
            hi_str = zeros(subsize(degs(ind)+1),1);
            Terms = cell2mat(children(h_i));
            for index = 1:size(Terms,2)
                monomial = Terms(index);
                c_degree = getTermDegree(monomial,z,varnum);
                c_index = findindex(a_f,c_degree);
                c_coeff = coeffs(monomial,z);
                c_coeff;
                hi_str(c_index) = c_coeff;
            end
            hm = max(abs(hi_str));
            hi_str = hi_str/hm;
            hs{ind} = hi_str;
        end

        msize = subsize(t+1);
        ysize = size(a_f,1);
        Mysize = subsize(t+1);

        equalities = cell(size(Zeqs,2),1);

        cvx_begin sdp quiet
        variable y(ysize,1)
        expression My(msize,msize)

        for ind1 = 1:msize
            v1 = a_f(ind1,:);
            for ind2 = 1:msize
                v2 = a_f(ind2,:);
                v = v1+v2;
                npos = findindex(a_f,v);
                My(ind1,ind2) = y(npos);
            end
        end
        minimize 1
        subject to
        My>=0;
        y(1)==1;

        % This is for M(hy)
        Mhy = cell(size(hs,1),1);
        for ind = 1: size(hs,1)
            trunc = t-dj(ind);
            msubsize = subsize(trunc+1);
            expression Msub(msubsize,msubsize)
            hc = hs{ind};
            hclen = size(hc,1);
            for pos1 = 1:msubsize
                for pos2 = 1:msubsize
                    dd1 = a_f(pos1,:);
                    dd2 = a_f(pos2,:);
                    alpha = dd1+dd2;
                    term = 0;
                    for index = 1:hclen
                        beta = a_f(index,:);
                        ydeg = alpha+beta;
                        yindex = findindex(a_f,ydeg);
                        term = term+y(yindex)*hc(index);
                    end
                    Msub(pos1,pos2) = term;
                end
            end
            Mhy{ind} = Msub;
            norm(Msub)<=1e-12;
        end
        cvx_end

        yseq = y;
        Mt = full(My);

        allranks = zeros(1,t+1);
        for ind = 0:t
            Ms = Mt(1:subsize(ind+1),1:subsize(ind+1));
            crank = svd_rank(Ms,threshold,decay);
            allranks(ind+1) = crank;
        end

        allranks;

        % find s candidates
        scand = [];
        for s = 2:t+1
            if allranks(s-1) == allranks(s)
                scand = [scand,s];
            end
        end

        if size(scand,2)>0
            for sind = 1:size(scand,2)
                s = scand(sind); % We first initialize s
                Bs = greedySieve(Mt,s,a_f,subsize,threshold,decay);
                Mb = getMb(Mt,Bs);
                Px = getPx(Mt,Bs,a_f);

                [Chi,Soln_cand,ComMtx,sum_err,avg_err] = getChi(Px,Mb);

                % Extracting a solution from chi. We verify z and if z is
                % accurate enough we plug it back into x

                indicators = ones(1,varnum);
                xind = zeros(1,varnum);
                terminate = 0;
                cerr = inf;
                zp = zeros(1,size(z,2)); %z previous
                zeval = zeros(1,size(z,2));
                % Get the most accurate soln.
                Bsize = size(Bs,2);
                % first, specify a candidate, then update
                csel = ones(varnum,1);
                
                while csel(1)<= Bsize
                    for ind = 1:varnum
                        zeval(ind) = Soln_cand(ind,csel(ind));
                    end
                    nerr = EvaluateErrors(hs,zeval,a_f);
                    if nerr<cerr
                        cerr = nerr;
                        zp = zeval;
                    end
                    csel = UpdateSelectedValue(csel,varnum,Bsize);
                end
                
                if cerr<= 1e-7 % This is NOT A GOOD THRESHOLD if we do not normalize G
                    flag = 0;
                    % Partial solution found, plug it back into x
                    for ind = 1:varnum
                        xsol(backindex(ind)) = zp(ind);
                    end
                    break;
                end
            end
        end

        td = td+1;

    end

end
%% Reconstruction

% Now, recover the solution:
n = size(Ap,1);
T = eye(n);
iT = eye(n);

p = 2;
q = 1;
Ac = Ap;
Bc = Bp;
Cc = Cp;
for ind = 1:n*(n-1)
    Ei = sym(eye(n));
    iEi = sym(eye(n));
    Ei(p,q) = xsol(ind);
    iEi(p,q) = -xsol(ind);
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

Ac = double(Ac);
Bc = double(Bc);
Cc = double(Cc);

% Round up the numbers in Ac, Bc and Cc. Set them to zero.
A = Ac;
B = Bc;
C = Cc;
for ind1 = 1:n
    for ind2 = 1:n
        if Asupp(ind1,ind2) == 0
            A(ind1,ind2) = 0;
        end
    end
end

for ind1 = 1:size(Bsupp,1)
    for ind2 = 1:size(Bsupp,2)
        if Bsupp(ind1,ind2) == 0
            B(ind1,ind2) = 0;
        end
    end
end

for ind1 = 1:size(Csupp,1)
    for ind2 = 1:size(Csupp,2)
        if Csupp(ind1,ind2) == 0
            C(ind1,ind2) = 0;
        end
    end
end

A
B
C

Aaug1 = [F,G*Cp;Bp*H,Ap];
Aaug2 = [F,G*C;B*H,A];


eC1 = eig(Aaug1);
eC2 = eig(Aaug2);






