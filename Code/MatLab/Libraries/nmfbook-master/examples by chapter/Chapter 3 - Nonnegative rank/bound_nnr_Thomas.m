% Applying the lower bounds presented in Chapter 3.3 on the slack matrix of
% the regular octagon 
clear all; clc; 
X = [1 1 0 0
     1 0 1 0
     0 1 0 1
     0 0 1 1]
[rc,geo,nnucnorm,tausos,hypsep] = lowerbounds_nnr(X,4); 