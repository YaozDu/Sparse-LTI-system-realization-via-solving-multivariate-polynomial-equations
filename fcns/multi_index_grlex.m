function [A,len] = multi_index_grlex(varnum, multdeg)
% Generate multi-indices in graded lex order:
% total degree increasing, and within each degree:
% higher powers of x1 come first
    len = [];
    A = [];

    for d = 0:multdeg
        A_d = gen_fixed_degree(varnum, d);
        A = [A; A_d];
        len = [len;size(A,1)];
    end
end
