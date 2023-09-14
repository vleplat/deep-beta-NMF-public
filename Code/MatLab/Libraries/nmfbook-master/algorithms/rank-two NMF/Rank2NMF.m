% Rank-two NMF computes a rank-two NMF of X
%
% If X >= 0 and rank(X) = 1, take W as any non-zero column of X and compute
% H accordingly 
% 
% If X >= 0 and rank(X) = 2, this is Algorithm 4.1 in Chapter 4, 
% with alpha1 = alpha2 = 0. This is an exact algorithm. 
% 
% Otherwise, the input matrix has is replaced by max(0,X2) where X2 is its 
% best rank-two approximation; cf. Chapter 6. 

function [W,H] = Rank2NMF(X) 

if min(X(:)) < 0 
    disp('The input matrix is not nonnegative: X <-- max(X,0).'); 
    X = max(0,X); 
end 
[u,s,v] = svds(X,3);
nX = norm(X,'fro'); 
% Case 0. rank(X) = 0
if s(1,1) <= eps
    disp('The input matrix is (close to) zero.');
    [m,n] = size(X); 
    W = zeros(m,1); 
    H = zeros(1,n); 
    return;
end
% Case 1. rank(X) = 1
if s(2,2) < 1e-9*s(1,1)
    disp('The input matrix has rank one, an optimal solution is returned.');
    sumX = sum(X); 
    [~,j] = max(sumX); 
    W = X(:,j); 
    H = sumX/sum(W); 
    return; 
end
% Case 2. rank(X) >= 2 
if s(3,3) > 1e-9*s(2,2)
	disp('The input matrix does not have rank two.');  
    disp('X <-- rank(X2,0) where X2 is the rank-2 truncated SVD of X.'); 
    disp('The solution computed might not be optimal.'); 
    [u,s,v] = svds(X,2);
    Xj = max(u*s*v',0); 
else
    disp('The input matrix has rank two, an optimal solution is returned.');
    Xj = X; 
end
% Remove zero columns
sumXj = sum(Xj); 
J = find(sumXj > 0); 
Xj = Xj(:,J);
% Normalize columns
[m,n] = size(Xj); 
for i = 1 : n
    Xj(:,i) = Xj(:,i)/sum(Xj(:,i));
end
% Compute extreme points of conv(Xj) hence a feasible solution W
% given that rank(X) = 2. 
W = zeros(m,2); 
[~,j1] = max( sum(Xj.^2) ); 
W(:,1) = Xj(:,j1); 
[~,j2] = max( sum( (Xj-repmat(Xj(:,j1),1,n)).^2) ); 
W(:,2) = Xj(:,j2); 
% Compute the corresponding weights
H = NNLS(W,X); 