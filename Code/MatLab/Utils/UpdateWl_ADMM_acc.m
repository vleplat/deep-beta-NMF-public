function [W,Z,U,residual,tcpu] = UpdateWl_ADMM_acc(W0,rho,maxIter,thres,innerloop,nu,epsi,rl,delta,alpha_l,lam_l,Wlminus1,Hl,Wp)
%UpdateWl_ADMM_Procedure Function:  ADMM-like procedure detailed in Section 3.2
% to tackle Problem (19) recalled here-under:

% min \lambda_\ell D(W_{\ell-1}, W H_\ell) + \lambda_{\ell+1} D(W, W_{\ell+1} H_{\ell+1}) + 
% \alpha_\ell \log\det\left( W^\top W + \delta I \right) 
% such that : W^\top e = e, W>=0
%% ------------------------------------
% Initialization
%--------------------------------------
% Init for Z and W
W = W0;
Z = W;
% Init dual variables
U = zeros(size(W0));
% Init of residual and loss fun arrays
residual = zeros(maxIter+1,1);
residual(1) = norm(W-Z,'fro');
verbose = 0;
mu = zeros(rl,1);
[m,n] = size(Wlminus1);
% Init extrapolation variables
Uchap = U;
Zchap = Z;
Z_prev = Z;
U_prev = U;
r = 10;

Y = inv(W'*W + delta*eye(rl));
Y_plus=max(0,Y);
Y_minus=max(0,-Y);
Eta=(8*alpha_l/lam_l*(W0*(Y_plus+Y_minus))+4*rho/lam_l*ones(m,rl)).*((Wlminus1./(W0*Hl+eps))*Hl'); 
Psi=4*alpha_l/lam_l*W0*(Y_plus+Y_minus)+2*rho/lam_l*ones(m,rl);

%% ----------------------------------
% Main optimization loop
%%-----------------------------------
tim = tic;
k = 0;
res = inf;
while k <= maxIter && res>thres
    %gamma_k update
    gammak = k/(k+r);
    % W-minimization
    for i=1:innerloop
     V = Zchap - Uchap;   
     % update mu (lagrangian multipliers)
     Phi=repmat(ones(1,n)*Hl',m,1)-4*alpha_l/lam_l*W0*Y_minus - rho/lam_l*V;
     mu=updatemu(Phi,Eta,Psi,W,mu,epsi);
     % update maxtrix W ("dictionaries")
     W = W0 .*(((Phi+ones(m,1)*mu').^2+Eta).^(1/2)-(Phi+ones(m,1)*mu'))./(Psi+eps);
     W = max(W,eps);
    end

    % Z-minimization
    V = W + Uchap;
    b=-log(Wp)-nu*V;
    Z = (lambertw(exp(-b)*nu))/(nu);

    % Dual updates
    resi = W - Z;
    U = Uchap + resi;
    
    % Extrapolation step
    Uchap = U + gammak*(U - U_prev);
    Zchap = Z + gammak*(Z - Z_prev);

    % Computation of ADMM residual
    residual(k+1) = norm(W-Z,'fro');
    res = residual(k+1);

    % Monitoring
    if k > 1 && mod(k,100)==0 && verbose==1
        fprintf('Iter: %d   | D_beta:  %0.2f  \n',k ,residual(k+1))
    end

    % Update previous iterates
    Z_prev = Z;
    U_prev = U;
    % Update iteration counter
    k = k + 1;


end
tcpu = toc(tim);

end