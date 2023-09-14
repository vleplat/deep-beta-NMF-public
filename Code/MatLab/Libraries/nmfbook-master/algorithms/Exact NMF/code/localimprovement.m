% From an initial solution, this code executes a NMF algorithm repeatedly
% until - a precision on the error is reached, or,
%       - the error is no longer decreasing enough, or,
%       - a time limit is exceeded.
%
% [W,H,e] = localimprovement(S,W,H,options)
%
% Input.
%    X             : (m x n) matrix to factorize
%   (W,H)          : initial matrices of dimensions (m x r) and (r x n)
%
%   The optional 'options' has the following fields:
%   algonmf        : see FroNMF.m, can be
%                       'HALS','MUUP','ANLS','ADMM','ALSH','FPGM'
%   tolerance      : the function stops if this precision on the error is obtained
%   convergiter    : number of iterations of an NMF algorithm between each test 
%   delta_t        : time of execution of an NMF algorithm between each test 
%   alpha          : if the error has not decreased by a factor alpha, the function stops
%   timemaxconverg : the function stops after 'timemaxconverg' seconds
% Output.
%   (W,H)          : nonnegative matrices s.t. WH approximates X
%    e             : final error ||M-WH||_F^2

function [W,H,e] = localimprovement(X,W,H,options)
    %addpath algosNMF/;
    if nargin <= 3
        algonmf        = 5;
        tolerance      = 1e-6;
        convergiter    = 1000;
        delta_t        = 1;
        alpha          = 0.99;
        timemaxconverg = 600;
    else       
        algonmf        = options.algonmf;
        tolerance      = options.tolerance;
        convergiter    = options.convergiter;
        delta_t        = options.delta_t;
        alpha          = options.alpha;
        timemaxconverg = options.timemaxconverg;
    end
    
    f0 = cputime;
    nX = norm(X,'fro');
    e0 = norm(X-W*H,'fro')/nX;
    while 1
        optionsNMFalgo.timemax = delta_t; 
        optionsNMFalgo.maxiter = convergiter; 
        init.W = W; init.H = H; 
        optionsNMFalgo.init = init; 
        optionsNMFalgo.algo = algonmf; 
        optionsNMFalgo.display = 0; 
        [W,H] = FroNMF(X,size(W,2),optionsNMFalgo); 
        %[W,H] = algoNMF(X,W,H,convergiter,delta_t,algonmf);
        e     = norm(X-W*H,'fro')/nX;
        if e<tolerance || e>e0*alpha || cputime-f0 > timemaxconverg
            break;
        end
        e0=e;
    end
end