% Tests for semi-NMF
clear all; clc; 
m = 20; 
n = 30; 
r = 5; 

% 1) The input matrix is nonnegative --> semi-NMF can be computed optimally
% with the SVD (gievn that XX^T is irreducible)
disp('**** Test semi-NMF on a nonnegative matrix ****'); 
X = rand(m,n); 
[W,H,e] = semiNMF(X,r); 
[u,s,v] = svds(X,r); 
fprintf('Error ||X-WH||_F of semi-NMF = %2.2f, error ||X-WH||_F of SVD = %2.2f.\n',...
                   e,norm(X-u*s*v','fro')); 
disp('press any button to continue'); 
pause; 

% 1) The input matrix is not nonnegative --> semi-NMF is harder to compute
disp('**** Test semi-NMF on a matrix with negative entries ****'); 
X = randn(m,n); 
[W,H,e] = semiNMF(X,r); 
[u,s,v] = svds(X,r); 
fprintf('Error ||X-WH||_F of semi-NMF = %2.2f, error ||X-WH||_F of SVD = %2.2f.\n',...
                   e(end),norm(X-u*s*v','fro')); 