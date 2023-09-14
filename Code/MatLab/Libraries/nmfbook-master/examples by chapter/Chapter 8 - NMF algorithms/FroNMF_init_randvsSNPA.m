% This code can be used to generate of Figure 8.5
% It compares two initializations: SNPA and random, to run A-HALS on the
% CBCL data set. 
clear all; close all; clc; 
% Data set
seed = 2020; 
load CBCL; 
r = 49; 
options.timemax = 5; 
% Initialization with SNPA (could use options.normalize = 1)
[K,W] = SNPA(X',r); 
H = X(K,:); 
W = W'; 
options.init.W = W; 
options.init.H = H;
% Paramters for A-HALS 
options.beta0 = 0; 
options.maxiter = Inf; 
options.timelimit = 100; 
options.algo = 'HALS';
options.alpha = 0.5; 
% Run A-HALS with SNPA vs. random init. 
disp('Running A-HALS with SNPA initialization...'); 
[Whe,Hhe,ehe,the] = FroNMF(X,r,options); 
[m,n] = size(X); 
options.init.W = rand(m,r); 
options.init.H = rand(r,n);
disp('Running A-HALS with random initialization...'); 
[Wher,Hher,eher,ther] = FroNMF(X,r,options); 
% Display
set(0, 'DefaultAxesFontSize', 25);
set(0, 'DefaultLineLineWidth', 2);
numarker = 20; 
semilogy(the,ehe,'-','MarkerSize', 15); hold on; 
semilogy(ther,eher,'--','MarkerSize', 15); hold on; 
grid on; 
h = legend('SNPA', 'rand');  %'Cuprite - SNPA', 'Cuprite - rand'); 
set(h,'Interpreter','latex'); 
ylabel('$\frac{\|X-WH\|_F}{\|X\|_F}$','Interpreter','latex'); 
xlabel('Time (s.)','Interpreter','latex'); 
axis([0 5 7.5*1e-2 0.5])