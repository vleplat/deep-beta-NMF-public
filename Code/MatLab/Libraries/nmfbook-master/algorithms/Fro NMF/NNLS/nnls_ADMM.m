function [H,WtW,WtX] = nnls_ADMM(X,W,options)

% Computes an approximate solution of the following nonnegative least
% squares problem (NNLS)
%
%           min_{H >= 0} ||X-WH||_F^2
%
% with ADMM, tackling the Lagrangian 
% 
% L(H,Y,Z) = ||X-WH||_F^2 + rho/2 ||H-Y||_F^2 + <Z, H-Y> 
% 
% via alternating optimization. 
%
% Input / Output; see nnls_input_output.m 
% 
% options.rho: parameter of the ADMM scheme; default: trace(WtW)/r

if nargin <= 2
    options = [];
end
if ~isfield(options,'delta')
    options.delta = 1e-6; % Stopping condition depending on evolution of the iterate V:
    % Stop if ||V^{k}-V^{k+1}||_F <= delta * ||V^{0}-V^{1}||_F
    % where V^{k} is the kth iterate.
end
if ~isfield(options,'inneriter')
    options.inneriter = 500; 
end

W = full(W); 
[m,n] = size(X);
[m,r] = size(W);
WtW = W'*W;
WtX = W'*X;

% If no initial matrices are provided, H is initialized as follows: 
if ~isfield(options,'init') || isempty(options.init)
    H = nnls_init(X,W,WtW,WtX); 
else
    H = options.init; 
end
% ADMM parameter 
if ~isfield(options,'rho')
    rho = trace(WtW)/r; 
else 
    rho = options.rho; 
end 
% Precompute the inverse 
invWtWrhoI = inv( WtW+ rho*eye(r) ); 
invWtWrhoIWtX = invWtWrhoI*WtX; 
% initialize Y and Z to zero 
Y = zeros(r,n); 
Z = zeros(r,n);
cnt = 1; 
Hp = Inf; 
while (norm(H-Y,'fro') > options.delta*norm(H,'fro') || norm(H-Hp,'fro') > options.delta*norm(H,'fro'))... 
            & cnt <= options.inneriter 
    Hp = H; 
    Y = max(0, H + Z/rho); 
    H = invWtWrhoIWtX + invWtWrhoI * (rho*Y - Z); 
    Z = Z + rho*(H - Y); 
    cnt = cnt + 1;
end
% Final projection to have V feasible 
H = max(H , 0);  
end % of function nnlsMU