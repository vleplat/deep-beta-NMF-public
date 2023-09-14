% UDISJ matrix 2^n x 2^n
clear all; clc; 

n = 3; 
if min(n) > 14
    warning('n is rather large... you might need to reconsider...'); 
end
A = dec2bin(0:2^n - 1) - '0'; 
B = A; 
X = zeros(2^n,2^n); 

for j = 1 : size(B,1) 
    X(:,j) = (1-A*B(j,:)').^2;
end

%% Compute an exact NMF with r = 7, it will fail
X 
r = 2^n - 1
options.maxiter     = 3;
[W,H] = exactNMFheur(X,r,options); 

% All lower bounds 
[rc,geo,nnucnorm,tausos,hypsep] = lowerbounds_nnr(X); 