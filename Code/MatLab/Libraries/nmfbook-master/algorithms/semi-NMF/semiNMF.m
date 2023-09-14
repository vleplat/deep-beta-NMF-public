function [W,H,e] = semiNMF(X,r,maxiter)
% [W,H,e] = BCDsemiNMF(X,r,maxiter)
%  This codes solves semi-NMF
% 
% min_{W, H} ||X-WH||_F such that H >= 0. 
%
% If the optimal rank-r approximation of X is semi-nonnegative, then the
% code returns an optimal solution. Otherwise, it uses block coordiante 
% descent for semi-NMF initialized with the SVD. 
%
% See N. Gillis and A. Kumar, "Exact and Heuristic Algorithms for Semi-
% Nonnegative Matrix Factorization", SIAM J. on Matrix Analysis and 
% Applications 36 (4), pp. 1404-1424, 2015. 
% See also https://sites.google.com/site/nicolasgillis/code for more
% details on this code and other optimization and initialization strategies. 
% 
%
% ****** Input ******
%   X      : m-by-n matrix 
%   maxiter: a maximum number of iterations
%
% ****** Output ******
%   (W,H)  : a semi-NMF of X \approx WH, H >= 0
%    e     : vector of the evolution of the error ||X-Xr||_F at each iteration k

% Step 1. truncated SVD and check whether it is semi-nonnegative
% This step is not recommended for sparse and large matrices because of the 
% explicit computation of Xr 
if prod(size(X)) < 1e8 && ~issparse(X)
    [u,s,v] = svds(X,r); 
    Xr = u*s*v'; 
    [seminnrank,W,H] = seminonnegativerank(Xr); 
else
    seminnrank = 0; 
end
% Main step
% if the truncated SVD is semi-nonnegative --> we can compute an optimal solution
if seminnrank == r % cf. Corollary 3.3 in the paper above.
    e = norm(X-Xr,'fro'); 
    return;
% Otherwise, use 2-BCD
else % SVD init + 2-BCD scheme
    [W,H] = SVDinitSemiNMF(X,r);
    if nargin <= 2
        maxiter = 100;
    end
    nX2 = norm(X,'fro')^2;
    for i = 1 : maxiter
        W = X/H; % Optimal solution
        options.init = H; 
        [H,WTW,WTX] = NNLS(W,X,options); % BCD on the rows of H
        if nargout >= 3
            e(i) = nX2 - 2*sum(sum( WTX.*H)) + sum(sum( (WTW).*(H*H') ) );
            e(i) = sqrt(max(0, e(i)));
        end
    end
end