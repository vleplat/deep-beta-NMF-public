% [W,H,e,t,etrue] = FroNMF(X,r,options) 
%
% This code allows you to run NMF algorithms based on the two-block 
% coordinate descent (2-BCD) method to solve the NMF problem
% 
%            min_{W >=0 and H >= 0} || X - WH ||_F^2 . 
% 
% It uses *extrapolation* on the standard 2-BCD scheme: 
% 0. Initialize (W,H)
% for k = 1, 2, ...
%   1. Update H
%   2. Update W
% end
% W and H are updated using an existing NNLS algorithm (see NNLSalgo.m),
% namely MU, HALS, active-set, FPGM, ADMM, ALS. 
% 
% Using options.beta0 = 0 allows you to run the algorithms without the 
% extrapolation. 
%
%
% *********
%   Input
% *********
% X : input m-by-n nonnegative matrix 
% r : factorization rank r 
%
% ---Options--- 
% .display : = 1, displays the evolution of the relative error (default), 
%            = 0 otherwise. 
% 
% .init    : init.W and init.H are the intial matrices. 
%           default: W = rand(m,r); H = rand(r,n)
% 
% .maxiter : maximum number of iterations. 
%           default: 500
% 
% .timemax : maximum allotated time
%           default: 60 seconds
%
% .accuracy: stops the algorithm if the relative error has not decreased by
%            an amount of options.accuracy durint 10 iterations (default = 1e-4)
%
% .algo    : = 'HALS', hierarchical alternating least suqares from Cichocki
%                     et al. (2008) -- default 
%            = 'ASET', the active set method from Kim and Park (2011)
%            = 'MUUP', the multiplicative updates from Lee and Seung (1999)
%            = 'FPGM', fast projected gradient method (Nesterov, 1983) 
%            = 'ADMM', alternating direction method of multipliers
%            = 'ALSH', alternating least squaures heuristic (projecting the
%                     unconstrained solution onto the nonnegative orthant
%% Parameters for the number if iterations of the NNLSalgo
% .alpha > 0       : maximum number of inner iterations depends on the
%                    cost between the first and the next updates: the 
%                    number of inner iterations is limited to 
%                    1 + alpha*this ratio. 
%                    default : 0.5
% .delta  in [0,1) : stops the inner iterations when the difference between
%                    two iterates is smaller than the difference of the
%                    first two iterates * delta. 
%                    defaut : 0.1 
% .inneriter >= 1  : maximum number of inner iterations for the NNLS 
%                    algorithm. 
%                    default : 100 
% 
%% Extrapolation parameters (we recommend to keep the default parameters)
% .extrapolprojH : = 1, extrapolation of H is performed after the update of W
%                  = 2, extrapolation of H is performed directly
%                  = 3, extrapolation of H is performed directly and it is
%                  projected onto the nonnegative orthant. (default) 
% .beta0         : initial value of the paramter beta un [0,1) 
%                  default: 0.5, except for MUUP where it is = 0. 
%         IF beta0=0 is used, the algorithms are run without extrapolation
% .eta           : decrease of beta when error increases 
%                  default: 1.5
% .gammabeta     : increase of beta when error decreases, it has to be < eta
%                  default: 1.1 for ANLS, 1.01 for A-HALS and MU
% .gammabetabar  : increase of the upper bound on beta when error decreases 
%                  It has to be < gammabeta. 
%                  default: 1.05 for ANLS, 1.005 for the others
%
% **********
%   Output
% **********
% (W,H) : rank-r NMF of X, that is, W is m-by-r, H is r-by-n, W>=0, H>=0, 
%          W*H approximates X. 
% e     : relative error of the iterates (W_y,H_n), that is,
%          ||X-W_yH_n||_F/||X||_F where W_y is the extrapolated sequence 
%          and H_n is the NNLS update of H_y. Note that W_y is not 
%          necessarily nonnegative (when options.extrapolprojH=2) but 
%          this value is cheap to compute. 
% t     : cputime --> plot(t,e) displays the evolution of the error e over
%           time
% etrue : if hp~=2, etrue = e; otherwise etrue is the error of (W_n,H_n)
%           which has to be recomputed from scratch hence the output etrue
%           should only be used for *benchmarking*. 
% 
% 
% Code modified from https://sites.google.com/site/nicolasgillis/code and 
% A.M.S. Ang and N. Gillis, Accelerating Nonnegative Matrix Factorization
% Algorithms using Extrapolation, Neural Computation 31 (2), pp. 417-439, 
% 2019. 


function [W,H,e,t,etrue] = FroNMF(X,r,options) 

cputime0 = cputime; 
[m,n] = size(X); 
if nargin < 3
    options = [];
end
if ~isfield(options,'display')
    options.display = 1; 
end
%% Parameters of NMF algorithm
if ~isfield(options,'init')
    W = rand(m,r); 
    H = rand(r,n); 
else
    W = options.init.W; 
    H = options.init.H; 
end
if ~isfield(options,'maxiter')
    options.maxiter = 500; 
end
if ~isfield(options,'algo')
    options.algo = 'HALS'; 
end
if ~isfield(options,'timemax')
    options.timemax = 60; 
end
if ~isfield(options,'accuracy')
    options.accuracy = 1e-4; 
end
% Parameters for the number of inner iterations: delta and alpha
if ~isfield(options,'delta')
    options.delta = 0.1; 
end
if ~isfield(options,'alpha')
    options.alpha = 0.5; 
end
if ~isfield(options,'inneriter')
    options.inneriter = 100; 
end
%% Parameters of extrapolation
if ~isfield(options,'extrapolprojH')
    options.extrapolprojH = 3; 
end
if ~isfield(options,'beta0')
    options.beta0 = 0.5; % options.beta0 = 0 is the standard NMF scheme
end
if options.algo == 'MUUP' % Extrapolation does not work well with MUUP
                          % since extrapolated variables might have
                          % negative entries
    options.beta0 = 0;
end
if ~isfield(options,'eta')
    options.eta = 1.5; 
end
if ~isfield(options,'gammabeta')
    if options.algo == 'ANLS'
        options.gammabeta = 1.1; 
    else
        options.gammabeta = 1.01; 
    end
end
if ~isfield(options,'gammabetabar')
    if options.algo == 'ANLS'
        options.gammabetabar = 1.05; 
    else
        options.gammabetabar = 1.005; 
    end
end
if options.eta < options.gammabeta || options.gammabeta < options.gammabetabar
    error('You should choose eta > gamma > gammabar.');
end
if options.beta0 > 1 || options.beta0 < 0
    error('beta0 must be in the interval [0,1].');
end 
H = full(H); 
W = full(W); 
%% Scale initialization so that argmin_a ||a * WH - X||_F = 1   
XHt = X*H'; 
HHt = H*H'; 
scaling = sum(sum(XHt.*W))/sum(sum( HHt.*(W'*W) )); 
W = W*scaling; 
%% Scale W and H so that columns/rows have the same norm, that is, 
% ||W(:,k)|| = ||H(k,:)|| for all k. 
normW = sqrt((sum(W.^2)))+1e-16;
normH = sqrt((sum(H'.^2)))+1e-16;
for k = 1 : r
    W(:,k) = W(:,k)/sqrt(normW(k))*sqrt(normH(k));
    d(k) = sqrt(normW(k))/sqrt(normH(k)); 
    H(k,:) = H(k,:)*d(k);
end
HHt = diag(d)*HHt*diag(d); 
XHt = XHt*diag(d);
%% Extrapolation variables 
Wy = W; 
Hy = H; 
%% Setting parameters and precomputations
beta(1) = options.beta0; 
betamax = 1; 
nX = norm(X,'fro'); 
e(1) = sqrt( nX^2 - 2*sum(sum( XHt.*W ) ) + sum( sum( HHt.*(W'*W) ) ) ) / nX; 
etrue(1) = e(1); 
emin = e(1); 
Wbest = W; Hbest = H; 
t(1) = cputime - cputime0; 
i = 1; 
if options.display == 1
    disp('Display of iteration number and relative error in percent:')
    if e(1) < 1e-4
        fprintf('%2.0f: %2.1d...', i, 100*e(1));
    else
        fprintf('%2.0f: %2.2f...', i, 100*e(1));
    end
    kdis = 2;
end
%% Main loop  
while i <= options.maxiter ... 
        && t(i) < options.timemax .... 
        && (i <= 12 || abs(e(i-1)-e(i-11))>=options.accuracy)
    %% NNLS for H 
    options.init = Hy; 
    [Hn,~,~] = NNLS(Wy, X, options); 
    % Extrapolation
    if options.extrapolprojH >= 2 
        Hy =  Hn + beta(i)*(Hn-H) ;
    else
        Hy = Hn;
    end
    % Projection
    if options.extrapolprojH == 3
        Hy = max(0, Hy); 
    end
    %% NNLS for W 
    options.init = Wy'; 
    [Wn,HyHyT,XHyT] = NNLS(Hy', X', options); 
    Wn = Wn'; XHyT = XHyT'; HyHyT = HyHyT'; 
    % Extrapolation
    Wy = Wn + beta(i)*(Wn-W);  
    if options.extrapolprojH == 1
        Hy =  Hn + beta(i)*(Hn-H) ;
    end
    %% Relative error 
    e(i+1) = max(0,nX^2 - 2*sum(sum( XHyT.*Wn ) ) + sum( sum( HyHyT.*(Wn'*Wn) ) ));
    e(i+1) = sqrt(e(i+1))/nX; 
    t(i+1) = cputime - cputime0; 
    %% Compute the 'true' error of the latest feasible solution (Wn,Hn)
    %% if options.extrapolprojH == 2, that is, if Hy is not feasible.
    %% This should be used (only) for benchmarking purposes.  
    if nargout >= 5
        nottakenintoaccount = cputime; 
        % etrue(i+1) is not equal to e(i) only when H_y is not feasible
        if options.extrapolprojH == 2 && min(Hy(:)) < 0 
        	etrue(i+1) = max(0,nX^2 - 2*sum(sum( (X*Hn').*Wn ) ) + sum( sum( (Hn*Hn').*(Wn'*Wn) ) ));
            etrue(i+1) = sqrt(etrue(i+1))/nX; 
        else 
            etrue(i+1) = e(i+1); 
        end
        % Do not take into account computation of etrue that are not needed
        t(i+1) = t(i+1) - max(0,cputime - nottakenintoaccount); 
    end
    %% Check the evolution of the error 
    if e(i+1) > e(i) 
        % If ADMM --> reduce delta
        if options.algo == 'ADMM'
            fprintf('\n');  
            options.delta = options.delta/10;
            options.inneriter = ceil(1.5*options.inneriter);  
            if options.display == 1
                disp('ADMM reduces the parameter delta and increase number of inner iterations.'); 
            end
        end
        % scale (W,H)
        normW = sqrt((sum(W.^2)))+1e-16;
        normH = sqrt((sum(H'.^2)))+1e-16;
        for k = 1 : r
            W(:,k) = W(:,k)/sqrt(normW(k))*sqrt(normH(k));
            H(k,:) = H(k,:)/sqrt(normH(k))*sqrt(normW(k));
        end
        % restart the scheme
        Wy = W;
        Hy = H;
        if i == 1
            betamax  = beta(i); 
        else
            betamax  = beta(i-1); 
        end
        beta(i+1) = beta(i)/options.eta; 
    else
        % keep previous iterate in memory
        W = Wn; 
        H = Hn; 
        beta(i+1) = min(betamax, beta(i)*options.gammabeta); 
        betamax = min(1, betamax*options.gammabetabar); 
    end
    %% Keep best iterate in memory
    if e(i+1) <= emin
        Wbest = Wn; 
        Hbest = max(Hy,0);
        emin = e(i+1); 
    end
    %% Display iteration number and error in percent
    if options.display == 1
        if i == 1
            % Display roughly every second, being a multiple of 10
            % with a lower bound of every 100 iterations
            displayparam = max(  1, 10^( log10( (t(1)+1e-2)^(-1)) )  ); 
            displayparam = max(  1, round(displayparam/10)*10 ); 
            displayparam = min(displayparam, 20); 
        end
        if i/displayparam == floor(i/displayparam) 
            if e(i+1) < 1e-4
                fprintf('%2.0f: %2.1d...', i, 100*e(i+1)); 
            else
                fprintf('%2.0f: %2.2f...', i, 100*e(i+1)); 
            end
            kdis = kdis + 1; 
        end
        if kdis == 11
            kdis = 1; 
            fprintf('\n'); 
        end
    end
    %% Increment
    i = i+1; 
end 
if options.display == 1 && mod(i,displayparam)
    if e(1) < 1e-4
        fprintf('%2.0f: %2.1d...', i, 100*e(i));
    else
        fprintf('%2.0f: %2.2f...', i, 100*e(i));
    end
end
if options.display == 1, fprintf('\n'); end
W = Wbest; 
H = Hbest; 
