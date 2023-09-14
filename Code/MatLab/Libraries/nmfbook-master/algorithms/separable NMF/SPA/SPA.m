function K = SPA(X,r,options)

% Successive projection algorithm for separable NMF
% 
% *** Description ***
% At each step of the algorithm, the column of R maximizing ||.||_2 is 
% extracted, and R is updated by projecting its columns onto the orthogonal 
% complement of the extracted column. The residual R is initializd with X. 
%
% See N. Gillis and S.A. Vavasis, Fast and Robust Recursive Algorithms 
% for Separable Nonnegative Matrix Factorization,  IEEE Trans. on Pattern 
% Analysis and Machine Intelligence 36 (4), pp. 698-714, 2014.
% 
% K = SPA(M,r,normalize) 
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
% .normalize : = 1, will scale the columns of X so that they sum to one,
%              hence matrix H will satisfy the assumption above for any
%              nonnegative separable matrix X. 
%              = 0, is the default value for which no scaling is
%              performed. For example, in hyperspectral imaging, this 
%              assumption is already satisfied and normalization is not
%              necessary. 
% .relerror   : stops the algorithm when 
%               max_k ||R(:,k)||_2 <= relerror * max_k ||X(:,k)||_2
%               where R is the residual, that is, R = X-X(:,K)H. 
%               default: 10^-12
% .display    : = 1, displays the iteration count (default)
%
% ****** Output ******
% K        : index set of the extracted columns. 
%
% This implementation of the algorithm is based on the formula 
% ||(I-uu^T)v||^2 = ||v||^2 - (u^T v)^2
% which allows to avoid the explicit computation of the residual
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
if ~isfield(options,'precision')
    options.precision = 1e-6; 
end
if options.normalize == 1
    % Normalization of the columns of M so that they sum to one
    D = spdiags((sum(X).^(-1))', 0, n, n); 
    X = X*D; 
end

normX0 = sum(X.^2); 
nXmax = max(normX0); 
normR = normX0; 

i = 1; 
% Perform r recursion steps (unless the relative approximation error is 
% smaller than 10^-9)
if options.display == 1
    fprintf('Extraction of the indices by SPA: \n'); 
end
while i <= r && sqrt(max(normR)/nXmax) > options.precision
    % Select the column of M with largest l2-norm
    [a,b] = max(normR); 
    % Norm of the columns of the input matrix M 
    % Check ties up to 1e-6 precision
    b = find((a-normR)/a <= 1e-6); 
    % In case of a tie, select column with largest norm of the input matrix M 
    if length(b) > 1, 
        [c,d] = max(normX0(b)); 
        b = b(d);
    end
    % Update the index set, and extracted column
    K(i) = b; 
    U(:,i) = X(:,b); 
    % Compute (I-u_{i-1}u_{i-1}^T)...(I-u_1u_1^T) U(:,i), that is, 
    % R^(i)(:,J(i)), where R^(i) is the ith residual (with R^(1) = M).
    for j = 1 : i-1
        U(:,i) = U(:,i) - U(:,j)*(U(:,j)'*U(:,i));
    end
    % Normalize U(:,i)
    U(:,i) = U(:,i)/norm(U(:,i)); 
    % Update the norm of the columns of M after orhogonal projection using
    % the formula ||r^(i)_k||^2 = ||r^(i-1)_k||^2 - ( U(:,i)^T m_k )^2 for all k. 
    normR = normR - (U(:,i)'*X).^2; 
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