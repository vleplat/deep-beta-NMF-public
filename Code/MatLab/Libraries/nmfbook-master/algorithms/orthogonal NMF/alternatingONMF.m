% [W,H,e] = alternatingONMF(X,r,options)
%
% This code allows you to solve ONMF based on the two-block 
% coordinate descent (2-BCD) method to solve the NMF problem
% 
%      min_{W and H} || X - WH ||_F^2 such that H >= 0 and HH^T = I_r.  
% 
% 0. Initialize (W,H)
% for k = 1, 2, ...
%   1. Update H
%   2. Update W
% end
% W and H are updated using closed-form solutions. 
% 
% Since the columns of W correspond to centroids of a subset of the columns
% of X, it will be nonnegative if X is. However, this code can also process
% input matrices, X, that are not nonnegative in which case W might not be
% nonnegative. 
% 
% *********
%   Input
% *********
% X : input m-by-n matrix 
% r : factorization rank r 
%
% ---Options--- 
% .display : = 1, displays the evolution of the relative error (default), 
%            = 0 otherwise. 
% 
% .init    : W = options.init is the intialization for W. 
%           default: W is computed via SNPA
% 
% .maxiter : maximum number of iterations. 
%           default: 100
%
% .delta   : stop algorithm when the relative error does not change by more
%            than options.delta between two iterations. 
%           default: 1e-4
%
% **********
%   Output
% **********
% (W,H) : rank-r ONMF of X, that is, W is m-by-r, H is r-by-n, H>=0, 
%          H^TH = I_r, W*H approximates X. 
% e     : relative error of the iterates (W,H), that is,
%          ||X-WH||_F/||X||_F 
% 
% This code is the 2-BCD method described in the paper 
% F. Pompili, N. Gillis, P.-A. Absil and F. Glineur, "Two Algorithms for 
% Orthogonal Nonnegative Matrix Factorization with Application to 
% Clustering", Neurocomputing 141, pp. 15-25, 2014. 

function [W,H,e] = alternatingONMF(X,r,options)

if nargin < 3
    options = [];
end
if ~isfield(options,'display')
    options.display = 1; 
end
if ~isfield(options,'maxiter')
    options.maxiter = 100;
end
if isfield(options,'init')
    W = options.init; 
else % Use SNPA 
    optionsSNPA.display = options.display; 
    if options.display == 1
        disp('Initialization by SNPA:'); 
    end
    K = SNPA(X,r,optionsSNPA); 
    W = X(:,K); 
    if length(K) < r
        error('SNPA was not able to extract r indices. This means that your data set does not even have r extreme rays.');
    end
end
if ~isfield(options,'delta')
    options.delta = 1e-6;
end 
[m,n] = size(X); 
[m,r] = size(W);
% Xn: normalized version of X, ||Xn(:,j)||_2  1 for all j
norm2x = sqrt(sum(X.^2,1)); 
Xn = X.*repmat(1./(norm2x+1e-16),m,1);  
normX2 = sum(X(:).^2); 
k = 1;
if options.display == 1
    disp('Iteration number and relative error of ONMF iterates:'); 
end
while k <= options.maxiter && (k <= 3 || abs(e(k-1)-e(k-2)) > options.delta)
    % H = argmin_H ||X-WH||_F, H >= 0, rows H orthogonal up to a scaling
    % of the rows of H
    H = orthNNLS(X,W,Xn); 
    % Normalize rows of H 
    norm2h = sqrt(sum(H'.^2,1))+1e-16;
    H = repmat(1./norm2h',1,n).*H;
    % W = argmin_W ||X-WH||_F, W >= 0
    W = X*H'; 
    % Compute relative error: 
    % ||X-WH||_F^2 = ||X||_F^2 - 2 <X,WH> + ||WH||_F^2
    %              = ||X||_F^2 - ||W||_F^2, since W = X*H^T & HH^T=I
    e(k) = sqrt( (normX2-sum(sum(W.^2)))/normX2 ); 
    if options.display == 1
        if e(k) < 1e-4
            fprintf('%2.0f: %2.1d...', k, 100*e(k));
        else
            fprintf('%2.0f: %2.2f...', k, 100*e(k));
        end
        if mod(k,10) == 0
            fprintf('\n');
        end
    end
    k = k + 1; 
end
if options.display == 1
    fprintf('\n'); 
end
end
