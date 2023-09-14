% Separable NMF applied on the CBCL data set; see Figure 7.2
clear all; clc; 

load CBCL; 
r = 49; 
options.normalize = 1; % Using this option is not necessary, but makes more 
                       % sense in the context of near-separable NMF because 
                       % there is no reason to believe that rows of W are 
                       % normalized (recall we apply SPA on X^T).   
K = SPA(X',r,options); % SNPA would work as well. 
H = X(K,:); 
W = NNLS(H',X')'; 
affichage(W,7,19,19);