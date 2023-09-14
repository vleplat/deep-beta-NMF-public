% Computing separable symmetric nonnegative matrix tri-factorization
% - Separable tri-symNMF - 
% 
% Given A = W S W^T where W^T is separable, this codes returns W and S 
% (up to permutation and scaling). 
% 
% Input  : symmetric matrix A (n-by-n) and r 
% Output : W (n-by-r) and S(r-by-r) such that WSW^T is close to A. 
% 
% This algorithm is implemented following the paper 
% Arora, Ge, Halpern, Mimno, Moitra, Sontag, Wu, Zhu, A practical algorithm 
% for topic modeling with provable guarantees, In International Conference 
% on Machine Learning (ICML), pp. 280-288, 2013. 

function [W,S] = septrisymNMF(A,r)

n = length(A); 
%% Idendity K such that A(K,K) = W(K,:) S W(K,:)^T 
options.normalize = 1; 
options.display = 0; 
K = SPA(A,r,options); 
%% Solve A(K,K) z = q = A(K,:)e 
q = A(:,K)'*ones(n,1); 
% y = A(K,K)\q; % This works in noiseless conditions
y = NNLS(A(K,K),q); % You might want to use other objective functions than 
                    % least squares depending on the noise statistic. 
%% Recover S and W 
S = diag(y)*A(K,K)*diag(y); 
% W = A(:,K)/( diag(z)*A(K,K) ); % This works in noiseless conditions
W = NNLS((diag(y)*A(K,K))', A(:,K)')'; 