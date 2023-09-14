% Application of WNMF for a toy recommender system example (6-by-5 matrix); 
% see Section 9.5 
% 
% Note that this problem is highly sensitive to initialization, even the
% rank-one problem is NP-hard (N. Gillis and F. Glineur, "Low-Rank Matrix 
% Approximation with Weights or Missing Data is NP-hard", SIAM J. on Matrix 
% Analysis and Applications 32 (4), pp. 1149-1165, 2011),  
% and hence you will not get the same solution as in the book. 
clear all; clc; 
% Zeros correspond to missing entries
X = [2 3 2 0 0;  
     0 1 0 3 2; 
     1 0 4 1 0; 
     5 4 0 3 2; 
     0 1 2 0 4; 
     1 0 3 4 3]
r = 3; 
P = X; 
P(P>0) = 1; 
options.r = r; 
rng(2020); 
options.nonneg = 1;    % Nonnegativity constraints in WLRA --> Weighted NMF
% We reinitialization WLRA as long as the solution is not within [0.5,6.5]
tet = 1; 
while tet >= 1 && tet <= 100
    [W,H,e] = WLRA(X,P,options); 
    WH = W*H'; 
    tet=tet+1; 
    if min(WH(:)) >= 0.5 && max(WH(:)) <= 6.5
        tet = 0;
    end
end
H = H'; 
% Normalization to have W in [0,5]
for k = 1 : r
    sck = max(W(:,k))/5; 
    W(:,k) = W(:,k)/sck; 
    H(k,:) = H(k,:)*sck;
end
W, H, 
% Relative residual error
fprintf('RMSE = ||X - W*H||_P/||P||_P = %2.2d%%.\n', ... 
    sqrt(sum(sum( ((X-W*H).^2).*P))) / sum( P(:) ) ); 
disp('Approximation:') 
W*H