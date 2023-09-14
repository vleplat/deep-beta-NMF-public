% This codes solves the weighted low-rank matrix approximation problem. 
% Given X (m x n), a nonnegative weight matrix P (m x n), and r or an 
% initial matrix W (mxr), it iteratively solves
%
% min_{W(mxr), H(nxr)}  ||X-WH^T||_P^2 
%                               + lambda (||W||_F^2+||H||_F^2), 
% 
% where ||X-WH^T||_P^2 = sum_{i,j} P(i,j) (X-WH^T)_{i,j}^2, 
% 
% using a block coordinate descent method (columns of W and H are the 
% blocks, like in HALS for NMF). 
% 
% It is possible to requires (W,H) >= 0 using options.nonneg = 1. 
% 
% ---Input---
%  X        : (m x n) matrix to factorize
%  P        : (m x n) nonnegative weight matrix
% *** options *** 
%  r        : factorization rank r (defaut = 1) 
%  nonneg   : if nonneg = 1: W>=0 and H>=0 (default=0)
%  W        : (m x r) initialization for factor W (default=rand)
% init      : if init == 1: greedy approach (default) 
%                init == 2: random init 
%                init == 3: SVD/NMF init (using 0's for the unknown entries) 
%  maxiter  : number of iterations (defaut=100) 
%  lambda   : penalization parameter (default=1e-6)
%
% ---Output---- 
% WH^T approximates X
% e(i) = relative error at iteration i =  ||X-WH^T||_P/||X||_P

function [W,H,e] = WLRA(X,P,options)

if P == 0
    error('P should contain at least one positive entry');
end

[m,n] = size(X); 
if nargin <= 2
    options = [];
end
if ~isfield(options,'r') && ~isfield(options,'W')
    options.r = 1;
end
if ~isfield(options,'W')
    options.W = rand(m,options.r); 
else
    options.r = size(options.W,2);
end
if ~isfield(options,'nonneg')
    options.nonneg = 0;
end
if ~isfield(options,'init') 
    options.init = 1;
end
if ~isfield(options,'lambda')
    lambda = 1e-6;
else
    lambda = options.lambda; 
end
if options.init == 1 % Greedy
    % Greedy initialization
    W = options.W; 
    R = X;
    for k = 1 : options.r
        Rp = P.*R;
        H(:,k) = (Rp'*W(:,k))./(P'*(W(:,k).^2)+lambda);
        if options.nonneg == 1
            H(:,k) = max(eps,H(:,k));
        end
        W(:,k) = (Rp*H(:,k))./(P*(H(:,k).^2)+lambda);
        if options.nonneg == 1
            W(:,k) = max(eps,W(:,k));
        end
        R = R - W(:,k)*H(:,k)';
    end
elseif options.init == 2
    W = options.W;
    H = rand(n,options.r);
    % Scaling
    normMnz = sqrt( sum(sum( (X.^2).*(P>0) ) ));
    normUVnz = sqrt( sum(sum( ((W*H').^2).*(P>0) ) ));
    H = 0.9*H/(normUVnz*normMnz+eps);
    % Scaling of rank-one factors
elseif options.init == 3
    if options.nonneg == 1
        [W,H] = FroNMF(X,options.r); 
        H = H';
    else
        [u,s,v] = svds(X,options.r); 
        W = u*sqrt(s); 
        H = v*sqrt(s); 
    end
end
[W,H] = scalingWH(W,H);
if ~isfield(options,'maxiter')
    options.maxiter = 100;
end  
error0 = sum(sum( (X.^2).*P ) ); 
R = X - W*H'; 
e(1) = sqrt(sum(sum( (R.^2).*P ) ) / (error0+eps)); 
% Main loop
for i = 2 : options.maxiter
    R = X - W*H'; 
    for k = 1 : options.r
        R = R+W(:,k)*H(:,k)'; 
        Rp = R.*P;
        W(:,k) = (Rp*H(:,k))./(P*(H(:,k).^2)+lambda);
        if options.nonneg == 1
            W(:,k) = max(eps,W(:,k)); 
        end
        H(:,k) = (Rp'*W(:,k))./(P'*(W(:,k).^2)+lambda);
        if options.nonneg == 1
            H(:,k) = max(eps,H(:,k)); 
        end
        R = R-W(:,k)*H(:,k)';
    end
    [W,H] = scalingWH(W,H); 
    e(i) = sqrt(sum(sum( (R.^2).*P ) ) / (error0+eps)); 
end

% Scaling of columns of W (mxr) and H (nxr) to have 
% ||W(:,k)|| = ||H(:,k)|| for all k
function [W,H] = scalingWH(W,H)
[m,r] = size(W); 
normW = sqrt((sum(W.^2)))+1e-16;
normH = sqrt((sum(H.^2)))+1e-16;
for k = 1 : r
    W(:,k) = W(:,k)/sqrt(normW(k))*sqrt(normH(k));
    H(:,k) = H(:,k)/sqrt(normH(k))*sqrt(normW(k));
end