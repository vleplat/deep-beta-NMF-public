% This code allows you to generate the Figures 7.5 and 7.6 that 
% compare the greedy near-separable NMF algorithms
clear all; clc;  

seed = 2020; 

m = 20;  % (in book: 40)
r = 10;  % (in book: 20)
n = 100; % (in book: 220)
condW = 6; 
nummat = 1; % Number of matrices generated per noise level (in book: 20)
% Choose tested algorithms where 
% 1=SPA,2=VCA,3=FastAnchorWords,4=SNPA,5=MVE-SPA,6=SPA-SPA
nalgo = 1:6; % You can select a subset of algorithms  
diri = 1; % Parameter for the Dirichlet distribution (in book: 0.5)
          % We use a larger value here because m is smaller. 
% Noise levels 
numnoiselevels = 10; % Number of noise levels (in book: 30)

% First experiment: well-conditioned Dirichlet
xp = 1; 
delta1 = linspace(0,1.3,numnoiselevels); 
rng(seed);  
[results1, timings1] = run_exp_greedy_sepNMF(m,n,r,xp,condW,delta1,nummat,nalgo,diri); 
scriptfigsepNMF(results1, timings1, delta1, nalgo, xp); 
% Second experiment: well-conditioned Middle-point
xp = 2; 
delta2 = linspace(0,0.7,numnoiselevels); 
rng(seed);  
[results2, timings2] = run_exp_greedy_sepNMF(m,n,r,xp,condW,delta2,nummat,nalgo,diri); 
scriptfigsepNMF(results2, timings2, delta2, nalgo, xp); 
% Third experiment: ill-conditioned Dirichlet
xp = 3; 
delta3 = [logspace(-7,1,numnoiselevels)]; 
rng(seed);
[results3, timings3] = run_exp_greedy_sepNMF(m,n,r,xp,condW,delta3,nummat,nalgo,diri); 
scriptfigsepNMF(results3, timings3, delta3, nalgo, xp); 
% Fourth experiment: ill-conditioned Middle-point
xp = 4;  
delta4 = [logspace(-3,log10(2),numnoiselevels)]; 
rng(seed);
[results4, timings4] = run_exp_greedy_sepNMF(m,n,r,xp,condW,delta4,nummat,nalgo,diri); 
scriptfigsepNMF(results4, timings4, delta4, nalgo, xp); 