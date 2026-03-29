function [empty_varspace,xdel] = detectIdealDim(Gb,x)
% x is in sorted order, so the smallest variable is the last index in x.
syms xtemp
xnum = max(size(x,1),size(x,2));
xmin = x(xnum); % So later when we reorder the terms, swap x whatever with this.
zerodegree = 1;
xdel = [];
for ind = 1:xnum
    % nGb = subs(Gb,x,y);
    nGb = subs(Gb,x(ind),xtemp);
    nGb = subs(nGb,xmin,x(ind));
    nGb = subs(nGb,xtemp,xmin);
    nGb = gbasis(nGb,x,'MonomialOrder','lexicographic');
    % Now, we use nGb to do variable detection:
    empty_varspace = 1;
    nGbsize = size(nGb,2);
    for ind2 = 1:nGbsize
        cvar = symvar(nGb(ind2));
        if isequal(cvar,xmin)
            empty_varspace = 0;
            break
        end
    end
    disp(x(ind));
    if empty_varspace == 1
        zerodegree = 0;
        xdel = [xdel,x(ind)];
        disp('Free');
    else
        disp('Not Free');
    end
    disp('---')
end
empty_varspace = zerodegree; %0-> positive dimension
end