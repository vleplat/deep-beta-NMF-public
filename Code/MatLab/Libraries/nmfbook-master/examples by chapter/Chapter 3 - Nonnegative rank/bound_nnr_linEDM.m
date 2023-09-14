% Comparing lower bounds for the linear EDMs of the form X(i,j) = (i-j)^2 

clear all; clc; 

n = 6; 
for i = 1 : n
    for j = 1 : n
        X(i,j) = (i-j)^2;
    end
end
X 

[rc,geo,nnucnorm,tausos,hypsep] = lowerbounds_nnr(X,n); 
% the restricted nonnegative rank of linear EDM is n 