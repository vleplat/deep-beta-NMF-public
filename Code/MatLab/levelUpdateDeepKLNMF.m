% Solve the problem
%       min_{W,H} KL(X,WH) + lam * KL(W,Wp)
% Using one step MU

function [W,H] = levelUpdateDeepKLNMF(H,X,W,Wp,lam,epsi,beta,HnormType)

[m,n] = size(X);
[~,r1] = size(W);
e = ones(m,n);

% Update of factor H using framework from (Leplat et al., 2021)
% % update mu (lagrangian multipliers)
prod = W*H;
JN1 = ones(n,1);
Jr1  = ones(r1,1);
Wt = W';
if beta == 1
    C = (Wt*(X./prod));
    D = Wt*e;
elseif beta==3/2
    C=(Wt*(((prod).^(beta-2)).*X));
    D=Wt*(prod).^(beta-1);
end
if strcmp(HnormType,'rows')
    [~,I] = min(D,[],2);
    idx = sub2ind(size(D),1:r1,I');
    mu_0_H = (D(idx)-C(idx).*(H(idx)))';
    mu_H = updatemu_hrows(C,D,H,1,mu_0_H,epsi);

    % % update matrix  H Coefficients ("activations")
    H = H .* (C./(D-mu_H*JN1'+eps));

elseif strcmp(HnormType,'cols')
    [~,I] = min(D,[],1);
    idx=sub2ind(size(D),I,1:n);
    mu_0_H = (D(idx)-C(idx).*(H(idx)))';
    mu_H = updatemu_hcols(C,D,H,1,mu_0_H,epsi);
    
    % % update matrix  H Coefficients ("activations")
    H = H .* (C./(D-Jr1*mu_H'+eps));
end
H = max(H,eps); 

% Update of factor W
Ht = H'; 
if beta == 1
    a = e*Ht - lam*log(Wp);
    b = W.*((X./(W*H)*Ht));
    W = max(eps,1/lam * b./(lambertw(1/lam*b.*exp(a/lam))));
elseif beta==3/2
    prod_pow = sqrt((W*H));
    W_pow = sqrt(W);
    A = (1./W_pow).*(prod_pow*Ht) + 2*lam;
    B = W_pow.*((X./prod_pow)*Ht);
    C = 2*lam*sqrt(Wp);
    W = 1/4*((C + sqrt(C.^2 + 4*A.*B))./(A)).^2;
    W = max(eps,W);
end
