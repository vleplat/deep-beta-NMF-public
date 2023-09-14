function [H,WtW,WtX] = nnls_HALSupdt(X,W,options)

% Computes an approximate solution of the following nonnegative least
% squares problem (NNLS)
%
%           min_{H >= 0} ||X-WH||_F^2
%
% with an exact block-coordinate descent scheme.
%
% See N. Gillis and F. Glineur, Accelerated Multiplicative Updates and
% Hierarchical ALS Algorithms for Nonnegative Matrix Factorization,
% Neural Computation 24 (4): 1085-1105, 2012.
%
% Input / Output; see nnls_input_output.m 
%
% Code modified from https://sites.google.com/site/nicolasgillis/code

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

eps0 = 0; cnt = 1; epsi = 1;
while epsi >= (options.delta)^2*eps0 && cnt <= options.inneriter %Maximum number of iterations
    nodelta = 0; 
    for k = 1 : r
        deltaH = max((WtX(k,:)-WtW(k,:)*H)/(WtW(k,k)+1e-16),-H(k,:));
        H(k,:) = H(k,:) + deltaH;
        % Last two line equivalent to the HALS update: 
        % H(k,:) = max(0, (WtX(k,:)-WtW(k,:)*H+WtW(k,k)*H(k,:))/(WtW(k,k)+1e-16) ) ; 
        nodelta = nodelta + deltaH*deltaH'; % used to compute norm(V0-V,'fro')^2;
        if H(k,:) == 0 % safety procedure 
            H(k,:) = 1e-16; 
        end 
    end
    if cnt == 1
        eps0 = nodelta;
    end
    epsi = nodelta;
    cnt = cnt + 1;
end

end % of function nnlsHALSupdt