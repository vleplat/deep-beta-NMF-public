function [H,WtW,WtX] = nnls_MU(X,W,options)

% Computes an approximate solution of the following nonnegative least
% squares problem (NNLS)
%
%           min_{H >= 0} ||X-WH||_F^2
%
% using multiplicative updates.  
%
% Input / Output; see nnls_input_output.m  
% 
% + options.epsilon gives lower bound on the entries of H for better
%   numerical performance; 
%     default = 2^(-52) (Matlab epsilon machine)
% 
% Remark. As opposed to the otehr NNLS algorithms, X, W and H need to be 
%         nonnegative. However, the MU can be adapted in this case by using
%         the decomposition X = X+ - X- and W = W+ - W- where X+, X-, W+
%         and W- are nonnegative, and working out the denominator and
%         numerator of the MU accordingly. 

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
if ~isfield(options,'epsilon')
    options.epsilon = 2^(-52);
end

W = full(W); 
[m,n] = size(X);
[m,r] = size(W);
WtW = W'*W;
WtX = W'*X;

if min(X(:)) < 0
    warning('The matrix X should be nonnegative. Zero entries set to 0.'); 
    X = max(X, 0); 
end
if min(W(:)) < 0
    warning('The matrix W should be nonnegative. Zero entries set to 0.'); 
    W = max(W, 0); 
end
if ~isfield(options,'init')  || isempty(options.init)
    H = nnls_init(X,W,WtW,WtX);
else
    H = full(options.init);
    if min(H(:)) == 0
        H(H==0) = 0.001*max(H(:)); % This is important for MU to be able to modify such entries
    end
end
if min(H(:)) < 0
    warning('The matrix H should be nonnegative. Zero entries set to 0.001.'); 
    H = max(H, 0.001); 
end
eps0 = 0; cnt = 1; eps = 1;
while eps >= options.delta*eps0 && cnt <= options.inneriter %Maximum number of iterations
    Hp = H;
    H = max(options.epsilon, H.*(WtX)./(WtW*H));
    if cnt == 1
        eps0 = norm(Hp-H,'fro');
    end
    eps = norm(Hp-H,'fro');
    cnt = cnt + 1;
end
end % of function nnlsMU