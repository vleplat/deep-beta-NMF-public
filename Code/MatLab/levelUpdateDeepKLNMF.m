% Solve the problem
%       min_{W,H} KL(X,WH) + lam * KL(W,Wp)
% Using one step MU

function [W,H] = levelUpdateDeepKLNMF(H,X,W,Wp,lam,epsi)

[m,n] = size(X);
[~,r1] = size(W);
e = ones(m,n);

% Update of factor H using framework from (Leplat et al., 2021)
% % update mu (lagrangian multipliers)
prod = W*H;
JN1 = ones(n,1);
Wt = W';
C = (Wt*(X./prod));
D = Wt*e;
[~,I] = min(D,[],2);
idx = sub2ind(size(D),1:r1,I');
mu_0_H = (D(idx)-C(idx).*(H(idx)))';
mu_H = updatemu_hrows(C,D,H,1,mu_0_H,epsi);

% % update matrix  H Coefficients ("activations")
H = H .* (C./(D-mu_H*JN1'+eps));
H = max(H,eps); 

% Update of factor W
Ht = H'; 
a = e*Ht - lam*log(Wp);
b = W.*((X./(W*H)*Ht));
W = max(eps,1/lam * b./(lambertw(1/lam*b.*exp(a/lam))));
