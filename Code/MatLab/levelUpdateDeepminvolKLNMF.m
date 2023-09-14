% Solve the problem
%       min_{W,H} \lambda_\ell * KL(X,WH) + \lambda_{\ell+1} * KL(W,Wp) + \alpha_\ell*log(WtW+delta*I)
%       s.t. et*W=et
% Using one step MU

function [W,H,res] = levelUpdateDeepminvolKLNMF(H,X,W,Wp,options,l)

% size computation
[m,n] = size(X);
[~,r] = size(W);

% loading parameters
lam_l = options.lambda(l);
epsi = options.epsi;
delta = options.delta(l);
alpha_l = options.alpha(l);
rho = options.rho;
thres = options.thres;
innerloop = options.innerloop;
maxIterADMM = options.maxIterADMM;

% Relative weight nu for Z-min. step
nu = rho/options.lambda(l+1);

% Update of W using ADMM procedure
if options.accADMM
    [W,~,~,res,~] = UpdateWl_ADMM_acc(W,rho,maxIterADMM,thres,innerloop,nu,epsi,r,delta,alpha_l,lam_l,X,H,Wp);
else
    [W,~,~,res,~] = UpdateWl_ADMM(W,rho,maxIterADMM,thres,innerloop,nu,epsi,r,delta,alpha_l,lam_l,X,H,Wp);
end

% Update H 
prod = W*H;
E = ones(m,n);
Wt = W';
H = (H .* (Wt*((E./prod).*X))./(Wt*E+eps));
H = max(H,eps);

