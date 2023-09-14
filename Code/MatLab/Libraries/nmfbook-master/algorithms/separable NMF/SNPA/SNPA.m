function [K,H] = SNPA(X,r,options) 

% Successive Nonnegative Projection Algorithm (variant with f(.) = ||.||^2)
%
% *** Description ***
% At each step of the algorithm, the column of X maximizing ||.||_2 is 
% extracted, and X is updated with the residual of the projection of its 
% columns onto the convex hull of the columns extracted so far. 
% 
% See N. Gillis, "Successive Nonnegative Projection Algorithm for Robust 
% Nonnegative Blind Source Separation", SIAM J. on Imaging Sciences 7 (2), 
% pp. 1420-1450, 2014.
%  
% [K,H] = SNPA(X,r,options) 
%
% ****** Input ******
% X          : an m-by-n matrix X.  
%              Ideally admitting a near-separable factorization, that is, 
%              X = WH + N, where 
%              conv([W, 0]) has r vertices,  
%              H = [I,H']P where I is the identity matrix, H'>= 0 and its 
%              columns sum to at most one, P is a permutation matrix, and
%              N is sufficiently small. 
% r          : number of columns to be extracted. 
% Options
% .maxitn    : number of inner iterations of the fast gradient method to
%               compute H 
%               default = 200; 
% .normalize : = 1, will scale the columns of X so that they sum to one,
%              hence matrix H will satisfy the assumption above for any
%              nonnegative separable matrix X. 
%              = 0, is the default value for which no scaling is
%              performed. For example, in hyperspectral imaging, this 
%              assumption is already satisfied and normalization is not
%              necessary. 
% .proj       : = 1, will use the constraints H^T e <= e, H >= 0. 
%               = 0, simply use the constraints H >= 0 (faster, but comes
%               with no theoretical guarantees at this point) 
% .relerror   : stops the algorithm when 
%               max_k ||R(:,k)||_2 <= relerror * max_k ||X(:,k)||_2
%               where R is the residual, that is, R = X-X(:,K)H. 
%               default: 10^-12
% .display    : = 1, displays the iteration count (default)
%
% ****** Output ******
% K        : index set of the extracted columns. 
% H        : optimal weights, that is, H = argmin_{Y >= 0} ||X-X(:,K)Y||_F
%
% Code modified from https://sites.google.com/site/nicolasgillis/code

[m,n] = size(X); 
% Options 
if nargin <= 2
    options = [];
end
if ~isfield(options,'normalize')
    options.normalize = 0; 
end
if ~isfield(options,'display')
    options.display = 1; 
end
if ~isfield(options,'maxitn')
    options.maxitn = 200; 
end
if ~isfield(options,'relerr')
    options.relerr = 1e-6; 
end
if ~isfield(options,'proj')
    options.proj = 0; 
end
if options.normalize == 1
    % Normalization of the columns of X so that they sum to one
    D = spdiags(((sum(X)+1e-16).^(-1))', 0, n, n); 
    X = X*D; 
end
% Initialization
normX0 = sum(X.^2); 
nXmax = max(normX0); 
i = 1; 
normR = normX0; 
XtUK = []; 
UKtUK = []; 
% Perform r recursion steps (unless the relative approximation error is 
% smaller than 10^-9)
if options.display == 1
    fprintf('Extraction of the indices by SNPA: \n'); 
end
% Main loop  
while i <= r & sqrt(max(normR)/nXmax) > options.relerr
    % Select the column of the residual R with largest l2-norm
    [a,b] = max(normR); 
    % Check ties up to 1e-6 precision
    b = find((a-normR)/a <= 1e-6); 
    % In case of a tie, select column with largest norm of the input matrix X 
    if length(b) > 1, 
        [c,d] = max(normX0(b)); 
        b = b(d); 
    end
    % Update the index set, and extracted column
    K(i) = b; 
    U(:,i) = X(:,b); 
    % Update MtUJ
    XtUK = [XtUK, X'*U(:,i)]; 
    % Update UJtUJ
    if i == 1
        UtUi = [];
    else
        UtUi = U(:,1:i-1)'*U(:,i); 
    end 
    UKtUK = [UKtUK, UtUi ; UtUi', U(:,i)'*U(:,i)]; 
    % Update residual 
    if i == 1
        % Fast gradient method for min_{y in Delta} ||M(:,i)-M(:,J)y||
        H = nnls_FPGM(X,X(:,K),options); 
    else
        H(:,K(i)) = 0; 
        h = zeros(1,n); h(K(i)) = 1; 
        H = [H; h]; 
        options.init = H; 
        H = nnls_FPGM(X,X(:,K),options); 
    end
    % Update the norm of the columns of the residual without computing it
    % explicitely. 
    if i == 1
        normR = normX0 - 2*( (XtUK').*H ) + ( H.*(UKtUK*H) );
    else
        normR = normX0 - 2*sum( (XtUK').*H ) + sum( H.*(UKtUK*H) );
    end
    if options.display == 1
        fprintf('%2.0f...', i); 
        if mod(i,10) == 0
            fprintf('\n');
        end
    end
    i = i + 1; 
end
if options.display == 1
    fprintf('\n'); 
end

end % of function SNPA