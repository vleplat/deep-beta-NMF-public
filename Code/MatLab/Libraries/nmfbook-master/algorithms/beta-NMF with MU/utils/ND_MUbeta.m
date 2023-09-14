% Numerator and denominator in the MU for the beta-divergence for the
% factor H, that is, H should be updated as follows: 
% 
% H <-- H. * (N./D) 
%
% and the error before the update of H is e. 
%
% Code modified from https://sites.google.com/site/nicolasgillis/code of
% the paper N. Gillis, L. T. K. Hien, V. Leplat and V. Y. F. Tan, 
% "Distributionally Robust and Multi-Objective Nonnegative Matrix 
% Factorization", January 2019, http://arxiv.org/abs/1901.10757

function [N,D,e] = ND_MUbeta(X,W,H,beta); 

epsilon = 2^(-52); 
if beta == 1 % Kullback-Leibler divergence
    if issparse(X)
        XdWH = blockrecursivecompwiseprodsparselowrank(X,W,H,@xdy); 
    else
        XdWH = X./(W*H+epsilon); 
    end 
    N = W'*XdWH; 
    D = repmat(sum(W)',1,size(X,2))+epsilon; 
    % error 
    Xnnz = X(X>0); 
    XdWHnnz = XdWH(X>0); 
    if nargout >= 3
        e = sum( Xnnz.* log( XdWHnnz+epsilon ) ) - sum( X(:) ) + sum( sum( D.*H ) ); 
    end
elseif beta == 2 % Frobenius norm
    N = W'*X; 
    WtW = W'*W; 
    D = WtW*H+epsilon; 
    if nargout >= 3
        e = 0.5*(norm(X,'fro')^2 - 2* sum(sum( N.*H ) ) + sum(sum( (WtW).*(H*H') ) ));
    end
else % Other beta divergences -- cannot handle sparse
    WH = W*H+epsilon; 
    N = W' * ((WH+epsilon).^(beta-2) .* X); 
    D = W' * ( (WH+epsilon).^(beta-1) ); 
    if nargout >= 3
        e = betadiv(X,WH,beta); 
    end
end