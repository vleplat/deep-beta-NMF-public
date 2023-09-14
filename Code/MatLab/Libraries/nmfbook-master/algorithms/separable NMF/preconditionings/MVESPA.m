% MVE SDP-based Preconditioning 
% 
% *** Description ***
% First, the semidefinite programming based preconditioning for 
% near-separable matrices is computed. Then, the successive projection 
% algorithm (SPA) is run on the preconditioned matrix (see FastSepNMF.m). 
%
% See Algorithm 2 in N. Gillis and S.A. Vavasis, "Semidefinite Programming 
% Based Preconditioning for More Robust Near-Separable Nonnegative Matrix 
% Factorization", SIAM J. on Optimization 25 (1), pp. 677-698, 2015.
%
% Code downloaded from https://sites.google.com/site/nicolasgillis/code
% 
% [K,Q,u] = MVESPA(M,r,options)
%
% ****** Input ******
% M = WH + N : a (normalized) noisy separable matrix, that is, W full rank, 
%              H = [I,H']P where I is the identity matrix, H'>= 0 and its 
%              columns sum to at most one, P is a permutation matrix, and
%              N is sufficiently small. 
% r          : number of columns to be extracted. 
% Options
% .normalize  : normalize=1 will scale the columns of M so that they sum to one,
%              hence matrix H will satisfy the assumption above for any
%              nonnegative separable matrix M. 
%              normalize=0 is the default value for which no scaling is
%              performed. For example, in hyperspectral imaging, this 
%              assumption is already satisfied and normalization is not
%              necessary. 
%              (default = 0.)
% .dimred     : parameter allowing to choose the dimensionality reduction
%              technique used (for m > r): 
%              dimred=0: use SPA to reduce dimensionality via pinv(X(:,K))*
%              where K is computed via SPA/ 
%              dimred=1: use the truncated SVD, that is, svds(M,r). 
%              dimred=2: use the truncated SVD, but computation based on MM'.
%              dimred~=1,2: use a random projection (u = randn(m,r))
%              (default = 0.)
% .display    : =1 displays the differents steps with running times  
%              affi=0 no display 
%              (default = 1.)
% 
% ****** Output ******
% K        : index set of the extracted columns. 
% Q        : the preconditioning (M <- Q*M).  
% u        : the linear dimensionality reduction (M <- u'*M). 

function [K,Q,u] = MVESPA(M,r,options)

[m,n] = size(M); 
if r > min(m,n)
    error('r has to be smaller than m and n'); 
end
% Options 
if nargin <= 2
    options = [];
end
if ~isfield(options,'normalize')
    options.normalize = 0; 
end
if ~isfield(options,'display')
    options.display = 1; 
end
if ~isfield(options,'dimred')
    options.dimred = 0; 
end

if options.normalize == 1
    % Normalization of the columns of M so that they sum to one
    D = spdiags((sum(M).^(-1))', 0, n, n); M = M*D; 
end

% Linear dimensionality rediction
etim = cputime; 
if options.display == 1
    fprintf('Linear dimensionality reduction started...\n')
end
[Mr, u] = lindimred(M,r,options.dimred); 
if options.display == 1
    fprintf(' Done. \n')
    fprintf('The linear dimensionality reduction took %2.2f seconds. \n',cputime-etim); 
end

% Preconditioning based on minimum volume ellipsoid 
Q = minvolell(Mr,1e-6,r,options.display); 

% Preconditioned SPA
optionsSPA.display = options.display; 
K = SPA(Q*Mr,r,optionsSPA); 