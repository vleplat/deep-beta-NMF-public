% semi-NMF on the CBCL data set; see Figure 5.4 

load CBCL; 
r = 49; 
% W nonnegative, not H: min_{W,H} ||X-WH||_F s.t. W >= 0 
[H,W,e] = semiNMF(X',r);   
W = W';
affichage(W,7,19,19); 
title('Basis elements of semi-NMF'); 