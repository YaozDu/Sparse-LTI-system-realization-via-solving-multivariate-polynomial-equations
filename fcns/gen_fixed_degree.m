function A = gen_fixed_degree(n, d)
% Generate all alpha in N^n with sum(alpha)=d
% in lex order: a1 descending

    if n == 1
        A = d;
        return;
    end

    A = [];

    for a1 = d:-1:0   % descending → gives your order
        B = gen_fixed_degree(n-1, d - a1);
        A = [A; [a1*ones(size(B,1),1), B]];
    end
end