% This codes allows you to generate the Figure 7.4 that compares  
% linear dimentionality reduction techniques used as a preprocessing for 
% MVE-SPA 
clear all; clc; 

seed = 2020; 
m = 20;  % in the book: 40
r = 10;  % in the book: 20
n = 100; % in the book: 200
condW = 6; 
diri = 0.5; 
nummat = 2;          % in the book: 10
numnoiselevels = 10; % in the book: 20 

nalgo = [5 5 5]; % This MVE-SPA x3 where the LDR is performed in different ways

delta = [logspace(-6,0,numnoiselevels)]; % in book: -5 instead of -6
xp = 3; 
[results, timings] = run_exp_greedy_sepNMF(m,n,r,xp,condW,delta,nummat,nalgo,diri); 
scriptfigsepNMF(results, timings, delta, nalgo, xp); 
legend('SPA','truncated SVD','random projection'); 