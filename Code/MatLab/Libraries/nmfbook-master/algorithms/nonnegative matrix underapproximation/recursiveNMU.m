% This code solves NMU, that is, 
% 
% min_{U>=0, V>=0} ||M-UV|| such that UV <= M, 
% 
% sequentially by computing one rank-one factor U(:,k)V(k,:) at a time. 
%
% see "Using Underapproximations for Sparse Nonnegative Matrix Factorization",
% N. Gillis and F. Glineur, Pattern Recognition 43 (4), pp. 1676-1687, 2010.
% and
% "Dimensionality Reduction, Classification, and Spectral Mixture Analysis 
% using Nonnegative Underapproximation", N. Gillis and R.J. Plemmons,
% Optical Engineering 50, 027001, February 2011.
% 
% website. http://sites.google.com/site/nicolasgillis/code
%
%
% [x,y] = recursiveNMU(M,r,Cnorm,maxiter)
%
% Input.
%   M              : (m x n) matrix to factorize.
%   r              : factorization rank, default = 1.
%   Cnorm          : Choice of the norm 1 or 2, default = 2.
%   maxiter        : number of iterations, default = 100.
%
% Output.
%   (U,V) : solution, UV^T "<=" M, U >= 0, V >= 0   ("." - approximately)

function [U,V] = recursiveNMU(M,r,Cnorm,maxiter)

[m,n] = size(M);
if nargin <= 3, maxiter = 100; end
if nargin <= 2, Cnorm = 2; end
if nargin <= 1, r = 1; end

% NMU Recursion
fprintf('Recursion started... \n')
for k = 1 : r
    % Initialization of (x,y) with an optimal rank-one NMF of M
    [u,s,v] = svds(M,1); x = abs(u)*sqrt(s); y = abs(v)*sqrt(s);
    U(:,k) = x; V(:,k) = y; 
    % Initialization of Lagrangian variable lambda
    R = M-x*y';
    lambda = max(0,-R);
    % Alternating optimization
    for j = 1 : maxiter
        % ***(x,y) Update***
        A = M-lambda;
        % l_1 norm minimization
        if Cnorm == 1
            x = max(0,(wmedian(A,y)));
            y = max(0,(wmedian(A',x)));
         % l_2 norm minimization   
         elseif Cnorm == 2 
            x = max(0,A*y); x = x/(max(x)+1e-16);
            y = max(0,(A'*x)/(x'*x));
        end
        % ***Lambda Update***
        if sum(x) ~= 0 && sum(y) ~= 0
            R = M-x*y'; U(:,k) = x; V(:,k) = y; 
            lambda = max(0,lambda - R/(j+1));
        else
            lambda = lambda/2;
            x = U(:,k); y = V(:,k); 
        end
    end
    M = max(0,M-x*y');
    if mod(k,10) == 0, fprintf('%1.0f...\n',k); 
    else fprintf('%1.0f...',k); end
    if k == r, fprintf('Done. \n',k); end
end


% WMEDIAN computes an optimal solution of
%
% min_x  || A - xy^T ||_1, y >= 0
%
% where A has dimension (m x n), x (m) and y (n),
% in O(mn log(n)) operations. Should be done in O(mn)...

function x = wmedian(A,y)

% Reduce the problem for positive entries of y
indi = y > 1e-16;
A = A(:,indi);
y = y(indi); 
[m,n] = size(A);
A = A./repmat(y',m,1);
y = y/sum(y);

% Sort rows of A, m*O(n log(n)) operations
[As,Inds] = sort(A,2);

% Construct matrix of ordered weigths
Y = y(Inds);

% Extract the median
actind = 1:m;
i = 1; 
sumY = zeros(m,1);
x = zeros(m,1);
while ~isempty(actind) % O(n) steps... * O(m) operations
    % sum the weitghs
    sumY(actind,:) = sumY(actind,:) + Y(actind,i);
    % check which weitgh >= 0
    supind = (sumY(actind,:) >= 0.5);
    % update corresponding x
    x(actind(supind)) = As(actind(supind),i);
    % only look reminding x to update
    actind = actind(~supind);
    i = i+1;
end