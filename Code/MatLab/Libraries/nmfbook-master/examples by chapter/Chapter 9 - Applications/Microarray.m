% Application of NMF on a gene microarray data set; see Section 9.4 
clear all; clc; 

load microarrayIFNbeta 
r = 3; 

options.maxiter = 100; 
[Ks,Hs] = SNPA(X,r); 
options.init.W = X(:,Ks); 
options.init.H = Hs; 
options.algo = 'HALS'; 
options.beta0 = 0; 

[W,H,e,t] = FroNMF(X,r,options); 

for i = 1 : r
    mwi =  max(W(:,i)); 
    W(:,i) = W(:,i)/mwi;
    H(i,:) = H(i,:)*mwi;
end

imagesc(W); 
colormap(gray); 
title('Basis matrix $W$', 'Interpreter', 'latex'); 