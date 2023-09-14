% This function aims to solve exactly or approximately the nonnegative
% least square (NNLS) problem of the form: 
% 
%         min_{H >= 0} ||X - WH||_F^2 
%
% using different algorithms: 
% (1) 'HALS': block coordinat descent method (on the rows of H)
% (2) 'ASET': active-set method from Kim and Park; 
%             see https://www.cc.gatech.edu/~hpark/nmfsoftware.html
% (3) 'MUUP': multiplicative updates (this is the approach
%       popularized by Lee and Seung but initially proposed in 
%       Daube-Witherspoon, M. E., & Muehllehner, G. (1986). An iterative 
%       image space reconstruction algorithm suitable for volume ECT. IEEE 
%       Trans. Med. Imaging, 5, 61–66. 
% (4) 'ADMM': alternating direction method of multipliers
% (5) 'FPGM': fast projected gradient method
% (6) 'ALSH': projection onto the nonnegative orthant of the unconstrained
%             least squares solution
% 
% *********
%   Input
% *********
% X : m-by-n matrix 
% W : m-by-r matrix 
%
% With MUUP, X and W need to be nonegative. 
%
% ---Options--- 
% .display : = 1, displays the evolution of the relative error (default), 
%            = 0 otherwise. 
% 
% .init      : intial iterate for H.  
%           default: init = [] meaning it will call nnls_init.m
% 
% .algo : = 'ASET', 'HALS', 'MUUP', 'FPGM', 'ADMM' or 'ALSH'; see above. 
% 
%% Parameters for the number if iterations of the NNLS algo
% .alpha > 0       : maximum number of inner iterations depends on the
%                    cost between the first and the next updates: the 
%                    number of inner iterations is limited to 
%                    1 + alpha*this ratio. 
%                    default : 0.5
% .delta  in [0,1) : stops the inner iterations when the difference between
%                    two iterates is smaller than the difference of the
%                    first two iterates * delta. 
%                    defaut : 1e-6
% .inneriter >= 1  : maximum number of inner iterations for the NNLS 
%                    algorithm. 
%                    default : 500 
%
% **********
%   Output
% **********
% H      :  approximate solution of min_{H >= 0} ||X - WH||_F^2
% WTW    : W'*W, and  
% WTX    : W'*X which can be used to compute the error efficiently since 
%           ||X - WH||_F^2 = ||X||_F^2 - 2 <WTX,H> + <WTW,H*H^T>
% 
% Code modified from https://sites.google.com/site/nicolasgillis/code

function [H,WTW,WTX] = NNLS(W,X,options) 

%% Options 
if nargin <= 2
    options = [];
end
if ~isfield(options,'algo')
    options.algo = 'HALS'; 
end
if ~isfield(options,'init')
    options.init= [];
end
%% number of innter iterations 
if ~isfield(options,'delta') 
    options.delta = 1e-6; 
end
if ~isfield(options,'inneriter') 
    options.inneriter = 500; 
end 
if isfield(options,'alpha') && sum(options.algo == 'ADMM') < 4
    [m,n] = size(X); 
    [m,r] = size(W); 
    alphainneriter = 1 + ceil( options.alpha*(nnz(X)+m*r)/(n*r) ); 
    options.inneriter = min(alphainneriter,  options.inneriter); 
end
%% Choose and run the algorithm
if options.algo == 'HALS' % Cyclic coordinate descent method; see for example 
    % Gillis, N., & Glineur, F. (2012). Accelerated multiplicative updates 
    % and hierarchical ALS algorithms for nonnegative matrix factorization. 
    % Neural computation, 24(4), 1085-1105.
    [H,WTW,WTX] = nnls_HALSupdt(X,W,options); 
elseif options.algo == 'ASET' % Alternating nonnegative least squares
    % Kim, J., & Park, H. (2011). Fast nonnegative matrix factorization: 
    % An active-set-like method and comparisons. 
    % SIAM Journal on Scientific Computing, 33(6), 3261-3281.
    % Sometimes, nnlsm_blockpivot fails and returns NaN solutions for rank
    % deficient matrices. This happens usually if a column of W is equal 
    % to zero, which we avoid by checking this case beforehand: 
    [H,~,WTW,WTX] = nnlsm_blockpivot( W, X, 0, max(0,options.init) ); 
    % If nnlsm_blockpivot fails for rank deficient matrices with
    % non-zero columns (this is extremely rare in practice), then we
    % execute nnlsm_activeset
    if sum( isnan(H(:)) ) > 0
        [H,~,WTW,WTX] = nnlsm_activeset( W, X, 0, 0, max(0,options.init) ); 
    end
elseif options.algo == 'FPGM' % Fast projected gradient method
    [H,WTW,WTX] = nnls_FPGM(X,W,options);  
elseif options.algo == 'MUUP' % Multiplicative updates
    [H,WTW,WTX] = nnls_MU(X,W,options); 
elseif options.algo == 'ADMM' % Alternating direction method of multipliers
    [H,WTW,WTX] = nnls_ADMM(X,W,options); 
    % Options will keep in memory the dual variables Z 
elseif options.algo == 'ALSH' % LS heuristic: project unconstrained solution
    WTW = W'*W; 
    WTX = W'*X; 
    r = size(W,2); 
    if cond(WTW) > 1e6
        delta = trace(WTW)/r; 
        H = max(0, solveNormalEqComb(WTW+1e-6*delta*eye(size(WTW)),WTX)); 
    else
        H = max(0, solveNormalEqComb(WTW,WTX)); 
    end
end