% Testing rank-two NMF code on a rank-two matrix, and a matrix close to
% having rank two 
clear all; clc; 
m = 10; 
n = 20; 
disp('----------------------'); 
disp('Case 1: X has rank two'); 
disp('----------------------'); 
X = rand(m,2)*rand(2,n); 
[W,H] = Rank2NMF(X); 
fprintf('The relative error ||X-WH||_F/||X||F is %2.2f%%.\n', ... 
    100*norm(X-W*H,'fro')/norm(X,'fro'))
disp('-------------------------------------'); 
disp('Case 2: X is close to have rank two'); 
disp('-------------------------------------'); 
X = rand(m,2)*rand(2,n) + 0.05*rand(m,n); 
[W,H] = Rank2NMF(X); 
fprintf('The relative error ||X-WH||_F/||X||F is %2.2f%%.\n', ... 
    100*norm(X-W*H,'fro')/norm(X,'fro'))