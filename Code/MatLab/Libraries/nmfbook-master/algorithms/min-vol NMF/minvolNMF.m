% Tihs codes solves minimim-volume NMF
% 
% min_{W,H} ||M-WH||_F^2 + lambda' * logdet(W^TW + delta I) 
% 
% where W >= 0, H >= 0, 
% and sum-to-one constraints on W or H: 
% H^T e <= e (model 1), or 
% H e    = e (model 2), or 
% W^T e  = e (model 3, default), or 
% H^T e  = e (model 4). 
% 
% This is solved by optimizing alternatively over W and H, 
% using a projected fast gradient method; see 
% V. Leplat, A.M.S. Ang, N. Gillis, "Minimum-Volume Rank-Deficient 
% Nonnegative Matrix Factorizations", ICASSP 2019, May 12-17, 2019, 
% Brighton, UK. 
% 
% ****** Input ****** 
% X : m-by-n matrix to factorize 
% r : factorization rank r 
% options: 
%   .model: = 1: H^T e  <=  e, cols of H sum to at most one (not wlog).  
%           = 2:   H e   =  e, rows of H sum to one (wlog). 
%           = 3:   W^T e =  e, columns of W sum to one (wlog) --> default 
%           = 4: H^T e   =  e, cols of H sum to one (not wlog) --> *new*  
%                There is a new numerical experiment using this last model;
%                see the folder "examples by chapter/Chapter 4 - Identifiability/minvolNMF_Moffet.m"
%                and also https://www.dropbox.com/s/iqp2saalnolfm2r/errataandmore_NMFbook.pdf?dl=0 
%   .lambda will be used to set up lambda' (default: 0.1), 
%       lambda' = lambda * ||M-WH||_F^2 / logdet( W^TW + delta I) 
%       where (W,H) is the initialization. 
%       This allows to balance the two terms in the objective (and also
%       make the algorithm unsensitive to scaling of the input matrix). 
%   .delta (default: 0.1) 
%   .maxiter (default: 100) 
%   .target: it allows to define a target for the relative error
%       ||X-WH||_F/||X||_F: lambda will be automatically tuned to attemp
%       achieving this target value (default: no target value). It is a 
%       simple heuristic: if ||X-WH||_F/||X||_F > target, lambda is 
%       decreased, and if ||X-WH||_F/||X||_F < target, lambda is increased. 
%   .(W,H): initialization (default: use of SNPA) 
%   .display: =1 displays the iteration count, = 0 no display
%               (default: 0).  
%
% ****** Outut ****** 
% (W,H) : low-rank approximation W*H of X, with H>=0, W>=0 of small volume,  
%         that is, logdet(W^TW + delta*I) is small, and  
%         if model = 1: H^T e  <=  e, or 
%                  = 2:   H e   =  e, or 
%                  = 3:   W^T e =  e (default), 
%                  = 4:   H^T e =  e. 
% e     : evolution of the error 
%         ||M-WH||_F^2 + lambda * logdet( W^TW + delta I) 
% err1  = evolution of ||M-WH||_F^2
% err2  = evolution of logdet( W^TW + delta I) 

function [W,H,e,err1,err2] = minvolNMF(X,r,options);

if nargin <= 2
    options = [];
end
if ~isfield(options,'model')
    options.model=3;
end
if ~isfield(options,'lambda')
    options.lambda=0.1;
end
if ~isfield(options,'delta')
    delta=0.1;
else
    delta=options.delta;
end
if ~isfield(options,'maxiter')
    options.maxiter=100;
end
if ~isfield(options,'inneriter')
    options.inneriter=10;
end
if isfield(options,'W') && isfield(options,'H') 
    W = options.W;
    H = options.H;
else
    if options.model == 1
        options.proj = 1;
    end
    [K,H] = SNPA(X,r,options); 
    W = X(:,K);
    optionsNNLS.init = H; 
    H = NNLS(W,X,optionsNNLS); 
    if length(K) < r
        warning('SNPA recovered less than r basis vectors.');
        warning('The data poins have less than r vertices.');
        r = length(K);
        fprintf('The new value of r is %2.0d.\n',r);
    end
end
if ~isfield(options,'display')
    options.display=1;
end 
% Normalization 
[W,H] = normalizeWH(W,H,options.model,X); 
% Initializations
normX2 = sum(X(:).^2);
normX = sqrt(normX2); 
WtW = W'*W;
WtX = W'*X;
% Initial error and set of of lambda
err1(1) = max(0,normX2-2*sum(sum(WtX.*H))+sum(sum( WtW.*(H*H'))));
err2(1) = log( det (WtW  + delta*eye(r) ) );  
lambda = options.lambda * max(1e-6,err1) / (abs( err2 ));
e(1) =  err1 + lambda * err2 ;
 % number of updates of W and H, before the other is updated
if options.display == 1
    disp('Iterations started: '); 
    fprintf('%1.0d ...', 1);
end
% projection model for H
if options.model == 1
    options.proj = 1;
elseif options.model == 2
    options.proj = 2;
elseif options.model == 3
    options.proj = 0;
elseif options.model == 4
    options.proj = 3;    
end
% Main loop 
for i = 2 : options.maxiter
    % *** Update W ***
    XHt = X*H';
    HHt = H*H';
    Y = inv( ( W'*W + delta*eye(r) ) );
    A = lambda*Y + HHt;
    if options.model <= 2 || options.model == 4
        W = FGMqpnonneg(A,XHt,W,options.inneriter,1); 
    elseif options.model == 3
        W = FGMqpnonneg(A,XHt,W,options.inneriter,2); 
    end
    % *** Update H ***
    options.init = H; 
    [H,WtW,WtX] = nnls_FPGM(X,W,options);
    err1(i) = max(0, normX2 - 2*sum(sum( WtX.*H ) )  + sum(sum( WtW.*(H*H') ) ) );
    err2(i) = log( det ( WtW + delta*eye(r) ) );
    e(i) = err1(i) + lambda * err2(i);
    if options.display == 1
        fprintf('%1.0d ...', i);
        if mod(i,10) == 0
            fprintf('\n');
        end
    end
    % Tuning lambda to obtain options.target relative error 
    if isfield(options,'target')
        if sqrt(err1(i))/normX > options.target+0.001
            lambda = lambda*0.95;
        elseif sqrt(err1(i))/normX < options.target-0.001
            lambda = lambda*1.05;
        end
    end
end