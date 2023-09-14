% Exact NMF of the nested hexagon problem 
clear all; clc; close all; 

%% Choose the parameter a 
a = 3; % a in (1, +Inf] 
%% Different cases depending on the value of a 
if a <= 1
    error('a should be chosen in (1, +Inf]'); 
elseif a < Inf
    Xa = [1     a     2*a-1  2*a-1  a  1
         1     1     a     2*a-1  2*a-1   a
         a     1     1     a  2*a-1   2*a-1
         2*a-1  a     1     1  a  2*a-1
         2*a-1  2*a-1  a     1  1  a
         a     2*a-1  2*a-1  a  1  1]  
    if a <= 2 % rank+(X) = 3
        r = 3; 
    elseif a <= 3 % rank+(X) = 4
        r = 4; 
    else % rank+(X) = 5 
        r = 5;
    end
else % a == Inf 
    Xa = [0 1 2 2 1 0
        0 0 1 2 2 1
        1 0 0 1 2 2
        2 1 0 0 1 2
        2 2 1 0 0 1
        1 2 2 1 0 0]
    r = 5; 
end
%% Compute an Exact NMF numerically
% exactNMFheur --> Heurtic algorithm for Exact NMF from
% A. Vandaele, N. Gillis, F. Glineur and D. Tuyttens, "Heuristics for
% Exact Nonnegative Matrix Factorization", Journal of Global Optimization
% 65 (2), pp 369-400, 2016.
% See the website https://sites.google.com/site/exactnmf/, 
% and the folder 'algorithms/Exact NMF' 
rng(1) 
options.heuristic   = 'ms2';
options.tolerance = 1e-9;
[W,H] = exactNMFheur(Xa,r,options) 
%% Display the solution in 2D, if possible 
% If rank(H)=3 and rank(W)=4, transpose solution to have a 2-D
% representation. This is possible because Xa is symmetric, up to
% permutation and scaling of its rows and columns. 
sw = svd(W); 
sh = svd(H); 
tol = 1e-6;  
if length(sw) > 3 && sw(4) > tol && sh(4) < tol
    Xa = Xa'; 
    Wold = W; 
    W = H'; 
    H = Wold'; 
    sw = sh; 
end
% Display NPP instance 
[P,U,V] = NPPrank3matrix(Xa); % display NPP instance corresponding to the 
                              % RE-NMF of a matrix Xa with rank(Xa)=3
% Display the solution W if it has rank 3 (if rank(W) > 3, conv(W) has 
% dimension strictly larger than 2)
if length(sw) == 3 || sw(4) < tol 
    sumWcol = sum(W); 
    W = W*diag(sumWcol.^(-1)); 
    H = diag(sumWcol)*H; 
    C = U\W;
    K = convhull(C(1,:),C(2,:)); 
    plot(C(1,K),C(2,K),'ks--','MarkerSize',15);
    legend('$\Delta^6 \cap $ col$(X_a)$', 'conv$(X_a)$', ... 
        'conv$(W)$', 'Interpreter','latex'); 
end