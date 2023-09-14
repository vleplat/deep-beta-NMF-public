% Given an m-by-r matrix W, comptues:
%
% gamma = min_i min_{x in Delta^r} ||W(:,i) - W(:,I) x||_2 / ||W(:,i)||_2
% where I = {1,2,...,r}\{i}
% 
% and k = index for which the minimum is achieved. 

function [gamma,k] = gamma2param(W) 

[m,r] = size(W); 
gamma = +Inf;
options.proj = 1; 
% Looping over the column of W
for i = 1 : r
    x = nnls_FPGM(W(:,i),W(:,[1:i-1 i+1:r]),options); 
    gammai =  norm(W(:,i) - W(:,[1:i-1 i+1:r])*x); 
    if gammai < gamma
        gamma = gammai;
        k = i;
    end
end