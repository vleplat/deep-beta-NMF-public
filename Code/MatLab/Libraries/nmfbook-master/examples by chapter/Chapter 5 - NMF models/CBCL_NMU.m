% sequential/recursive NMU on CBCL; see Figure 5.6 
clear all; clc; 
load CBCL; 
X = X';
[m,n] = size(X); 
r = 49; 
[W,H] = recursiveNMU(X,r,2,200); 
affichage(H,7,19,19); 