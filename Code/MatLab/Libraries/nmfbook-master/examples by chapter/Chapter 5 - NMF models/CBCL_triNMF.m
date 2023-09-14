% Performing tri NMF *sequentially* on the CBCL data set (this is a naive 
% approach): First compute X = WH, then W = W1*H1. 
% The decomposition W = W1*H1 is done with sparse NMF because otherwise the
% trivial decomposition is possible (since W1 has more columns than W); 
% see Figure 5.10 
clear all; clc; 
load CBCL 
r = 25; 
r1 = 49; 
options.maxiter = 500; 
options.timemax = Inf; 
%% First NMF: WH 
disp('*** First layer computation:  X = W*H - standard NMF  ***'); 
[W,H,e] = FroNMF(X,r,options);
%% Scale W and H so that columns/rows have the same norm, that is, 
% ||W(:,k)|| = ||H(k,:)|| for all k, 
% otherwise the second layer computation might be highly biased! 
% (Note that other scalings are possible, e.g., ||W(:,k)||=1 for all k.)
normW = sqrt((sum(W.^2)))+1e-16;
normH = sqrt((sum(H'.^2)))+1e-16;
for k = 1 : r
    W(:,k) = W(:,k)/sqrt(normW(k))*sqrt(normH(k));
    H(k,:) = H(k,:)/sqrt(normH(k))*sqrt(normW(k));
end
%% Second NMF: W = W1 H1
% We use sparse NMF with sparsity 0.85
% We need to use sparisty otherwise the factorization is highly ill-posed
options.sW = 0.85; % W1 is required to have 85% sparsity (Hoyer measure)
options.colproj = 1; 
disp('*** Second layer computation: W = W1*H1 - sparse NMF ***');
[W1,H1,e2] = sparseNMF(W,r1,options);
% Display results 
affichage(W1,7,19,19); title('$W_1$ in the triNMF $W_1 W_2 H \approx X$','Interpreter','latex');
affichage(W,5,19,19); title('$W_1 W_2$ in the triNMF $W_1 W_2 H \approx X$','Interpreter','latex');