% NMF with the Frobenius norm applied on the CBCL facial image data set 
clear all; clc;  
load CBCL; 
r = 49; 
rng(1234)
[W,H,e] = FroNMF(X,r); % see folder 'algorithms/Fro NMF' 
affichage(W,7,19,19);  % Display the basis images (columns of W) 
title('Basis elements: columns of $W$ reshaped as images', 'Interpreter', 'latex'); 