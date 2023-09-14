% Generates the slack matrix of the regular n-gon (hence the matrix
% has dimension n x n). 
% See `Heuristics for Exact Nonnegative Matrix Factorization', A. Vandaele,  
% N. Gillis, F. Glineur and D. Tuyttens, for more details. 

function [S,r]=slackregularngon(n)
    line1 = cn(n,0:n-1); 
    S     = zeros(n,n);
    shift = 1:n;
    for i=1:n
        S(i,shift) = line1;
        shift      = [shift(2:end) shift(1)];
    end
    if nargout>1
        r=nnrank(n);
    end
    S = max(S,0); 
end

function r=nnrank(n)
    k0=ceil(log2(n));
    k1=k0-1;
    k2=k0-2;
    
    lb1 = power(2,k1);
    ub1 = power(2,k1) + power(2,k2);
    lb2 = ub1;
    ub2 = power(2,k0);
    
    if lb1 < n && n <= ub1
        r = 2*k0 - 1;
    elseif lb2 < n && n <= ub2
        r = 2*k0;
    end
end

function c=cn(n,k)
    c = cos(pi/n)-cos(pi/n + 2*pi*k/n);
end

function d=dn(n,k)
    d=cn(n,k)/cn(n,1);
end