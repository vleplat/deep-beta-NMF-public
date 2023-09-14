function [k,normM,U] = SPArefine(M,K)

% Successive projection algorithm for separable NMF
% 
% Given K, what is the next best index? 
% 
% 
% *** Description ***
% At each step of the algorithm, the column of M maximizing ||.||_2 is 
% extracted, and M is updated by projecting its columns onto the orthogonal 
% complement of the extracted column. 
%
% See N. Gillis and S.A. Vavasis, Fast and Robust Recursive Algorithms 
% for Separable Nonnegative Matrix Factorization,  IEEE Trans. on Pattern 
% Analysis and Machine Intelligence 36 (4), pp. 698-714, 2014.
% 
% [J,normM,U] = SPArefine(M,r,normalize) 
%
% ****** Input ******
% M = WH + N : a (normalized) noisy separable matrix, that is, W is full rank, 
%              H = [I,H']P where I is the identity matrix, H'>= 0 and its 
%              columns sum to at most one, P is a permutation matrix, and
%              N is sufficiently small. 
% r          : number of columns to be extracted. 
% normalize  : normalize=1 will scale the columns of M so that they sum to one,
%              hence matrix H will satisfy the assumption above for any
%              nonnegative separable matrix M. 
%              normalize=0 is the default value for which no scaling is
%              performed. For example, in hyperspectral imaging, this 
%              assumption is already satisfied and normalization is not
%              necessary. 
%
% ****** Output ******
% J        : index set of the extracted columns. 
% normM    : the l2-norm of the columns of the last residual matrix. 
% U        : normalized extracted columns of the residual. 
%
% --> normM and U can be used to continue the recursion later on without 
%     recomputing everything from scratch. 
%
% This implementation of the algorithm is based on the formula 
% ||(I-uu^T)v||^2 = ||v||^2 - (u^T v)^2. 

[m,n] = size(M); 

normM = sum(M.^2); 
normM1 = normM; 

% Compute residual norms using successive projections 
for i = 1 : length(K); 
    U(:,i) = M(:,K(i)); 
    % Compute (I-u_{i-1}u_{i-1}^T)...(I-u_1u_1^T) U(:,i), that is, 
    % R^(i)(:,J(i)), where R^(i) is the ith residual (with R^(1) = M).
    for j = 1 : i-1
        U(:,i) = U(:,i) - U(:,j)*(U(:,j)'*U(:,i));
    end
    % Normalize U(:,i)
    U(:,i) = U(:,i)/norm(U(:,i)); 
    % Update the norm of the columns of M after orhogonal projection using
    % the formula ||r^(i)_k||^2 = ||r^(i-1)_k||^2 - ( U(:,i)^T m_k )^2 for all k. 
    normM = normM - (U(:,i)'*M).^2; 
    i = i + 1; 
end

% Index k to add to K 
[a,b] = max(normM); 
% Check ties up to 1e-6 precision
b = find((a-normM)/a <= 1e-6); 
% In case of a tie, select column with largest norm of the input matrix M 
if length(b) > 1, 
    [c,d] = max(normM1(b)); 
    b = b(d); 
end
k = b; 