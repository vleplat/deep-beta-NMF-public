% This codes allows to solve projective NMF
% 
%      min_{W >= 0} ||X - WW^TX||_F^2 
% 
% using multiplicative updates; see 
% Yang, Z., & Oja, E. (2010). Linear and nonlinear projective nonnegative 
% matrix factorization. IEEE Transactions on Neural Networks, 21(5), 
% 734-749. 
%
% Note: in the book, projective NMF is presented on the matrix transpose,
% that is, the model min_{H >= 0} ||X - XH^TH||_F^2 is presented which is 
% equivalent up to transposition of X. 
% 
% *** Input ***
% X        : data matrix 
% r        : factorization rank
% --- options ---
% .init    : initial matrix W 
% .maxiter : number of iterations
% .display : display the evolution of iterations
%
% *** Output ***
% W        : nonnegative matrix s.t. WW^TX is close to X
% e        : evolution of the relative error, ||X - WW^TX||_F/||X||_F

function [W,e] = projectiveNMF(X,r,options); 

[m,n] = size(X); 
% X: m by n 
% W: m by r
if nargin <= 2
    options = [];
end
if ~isfield(options,'init')
    W = rand(m,r);
else
    W = options.init; 
end
if ~isfield(options,'maxiter')
    options.maxiter = 200;
end
if ~isfield(options,'delta')
    options.delta = 1e-4;
end
if ~isfield(options,'display')
    options.display = 1;
end

if options.display == 1
    disp('Running projective NMF:'); 
end
nX2 = sum(sum(X.^2));
i = 1; 
Wp = 0; 
while i <= options.maxiter && norm(W-Wp,'fro') > options.delta*norm(W,'fro')
    Wp = W; 
    XtW = X'*W;   % n by r 
    XXtW = X*XtW; % m by r 
    WtW = W'*W;   % r by r
    XtWtXtW = XtW'*XtW; % r by r 
    % Optimal scaling of the initial solution
    % This is important otherwise the algorithm oscillates
    alpha = sum(sum(XXtW.*W)) / sum(sum(WtW.*XtWtXtW));
    W = W*sqrt(alpha);
    % Update the other factors accordingly
    XtW = sqrt(alpha)*XtW;
    XXtW = sqrt(alpha)*XXtW;
    WtW = alpha*WtW;
    XtWtXtW = alpha*XtWtXtW;
    if nargout >= 2
        e(i) = sqrt(max(0, nX2 - 2*sum(sum(XXtW.*W)) + sum(sum(WtW.*XtWtXtW)))); 
    end
    % MU by Yang and Oja
    W = W .* ( 2*XXtW ) ./ ( W*(XtWtXtW) + XXtW*(WtW) ); 
    % The order in which the multiplications are performed are chosen so as
    % to minimize the computational cost. 
    if options.display == 1
        fprintf('%2.0f...',i);
        if mod(i,10) == 0
            fprintf('\n');
        end
    end
    i = i + 1; 
end
e = e/norm(X,'fro'); 