%   This code solves the betaNMF problem
%
%          min_{W,H >= options.epsilon} D_{options.beta}(X,WH) 
%
%   where options.beta is the beta divergence used to measure the error, 
%   and options.epsilon is a small nonnegative constant (we recommend to
%   use the machine epsilon, which is the default value; cf. the discussion
%   about the zero locking phenomenon of MU in the book). 
% 
% 
% ****** Input ******
% X     :  the nonnegative input matrix pair 
% r     :  the rank of the sought approximation 
% 
% ---Options---    
% .maxiter    : the maximum number of iterations performed by the algorithm 
%             -default = 500. 
% .timemax   : the maximum time in seconds alloted to the algorithm 
%             -default = 60. 
% .beta       : beta divergence considered
% .epsilon    : lower bound on the entries of W and H to ensure convergence
%             -default: Matlab machine precision 2^-52
% .accuracy   : stop iterations if the relative error does not decrease by 
%               this value between two iterations (default: 1e-6) 
% .W and .H    : initial values for W and H. 
%                default: rand(m,r) and rand(r,n)
% .display in {0,1} : =1 displays the evolution of the objective (default), 
%                     =0 otherwise. 
% 
% ****** Output ******
% (W,H)     : W>=0 and H>=0, and WH approximates X according the the
%               criterion described above. 
% e         : e gives the values of the objective functions during the 
%             iterative process, that is, D_beta(X,WH)
%
% Code modified from https://sites.google.com/site/nicolasgillis/code of
% the paper N. Gillis, L. T. K. Hien, V. Leplat and V. Y. F. Tan, 
% "Distributionally Robust and Multi-Objective Nonnegative Matrix 
% Factorization", January 2019, http://arxiv.org/abs/1901.10757


function [W,H,e] = betaNMF(X,r,options);

[m,n] = size(X); 
if nargin <= 2
    options = [];
end
if ~isfield(options,'maxiter')
    options.maxiter = 500;
end
if ~isfield(options,'timelimit')
    options.timemax = Inf;
end
if ~isfield(options,'beta')
    warning('Beta not specified, default selected: beta=1 => KL divergence.'); 
    options.beta = 1;
end
if ~isfield(options,'epsilon')
    options.epsilon = 2^(-52);
end
if ~isfield(options,'accuracy')
    options.accuracy = 1e-4;
end
if ~isfield(options,'W')
    W = rand(m,r);
else
    W = max(options.epsilon,options.W); 
end
if ~isfield(options,'H')
    H = rand(r,n);
else
    H = max(options.epsilon,options.H); 
end
if ~isfield(options,'display')
    options.display = 1;
end

%if options.beta == 2
%    warning('Since beta=2, you might want to use more efficient algorithms; see FroNMF.m'); 
%end

i = 1; 
cpuinit = cputime; 
while i <= options.maxiter ... 
        && cputime <= cpuinit+options.timemax... 
        && (i <= 12 || abs(e(i-1)-e(i-11)) > options.accuracy*abs(e(i-2))) 
    % Update of (W,H) with the MU
    H = MUbeta(X,W,H,options.beta,options.epsilon); 
    [W,ei] = MUbeta(X',H',W',options.beta,options.epsilon); 
    W = W'; 
    % Error
    e(i) = ei; 
    % Scaling: the maximum entry in each column of W is 1
    for k = 1 : r
        mxk = max( W(:,k) ); 
        W(:,k) = W(:,k)/mxk;
        H(k,:) = H(k,:)*mxk;
    end 
    % Display evolution of the iterations
    if options.display == 1
        if mod(i,10) == 0
            fprintf('%1.2d...',i);
        end
        if mod(i,100) == 0
            fprintf('\n');
        end
    end
    i = i+1; 
end
if options.display == 1
    fprintf('\n');
end
