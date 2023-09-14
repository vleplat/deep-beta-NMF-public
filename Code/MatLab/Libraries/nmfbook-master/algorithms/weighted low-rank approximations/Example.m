% Example on a randomly generated matrix with about 50% missing entries
clear all; clc; close all; 
m = 100; 
n = 100; 
r = 5; 
W0 = randn(m,r); 
H0 = randn(n,r); 
X = W0*H0'; 
P = round(rand(m,n)); % ~50% missing entries, P is binary 

options.r = r; 
[W,H,e] = WLRA(X,P,options); 
% Relative residual error
fprintf('Relative error ||X-WH''||_F/||X||_F = %2.2f%%.\n', ... 
    100*norm(X - W*H','fro')/norm(X,'fro')); 
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultLineLineWidth', 2);
figure; 
semilogy(e); 
xlabel('Iterations'); 
ylabel('||X-WH''||_P / ||X||_P'); 