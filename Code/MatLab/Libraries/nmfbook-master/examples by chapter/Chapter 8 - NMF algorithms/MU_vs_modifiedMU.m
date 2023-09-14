% Compare the original MU with the modified MU that uses a positive lower 
% bound (namely, machine epsilon) for the entries of W and H; see 
% Figure 8.2. 
clear all; clc; 
% Generate sparse data set randomly 
m = 500;
n = 1000; 
r = 40; 
X = sprand(m,n,0.01); 
% Random initialization 
W = rand(m,r); 
H = rand(r,n); 
% Parameters of MU 
options.W = W; 
options.H = H; 
options.maxiter = 100; 
options.beta = 2; 
options.accuracy = 0; 
% Standard MU 
options.epsilon = 0; 
disp('Running MU'); 
[Ws,Hs,es] = betaNMF(X,r,options);
% Modified MU 
options.epsilon = eps; 
disp('Running modified MU'); 
[Wm,Hm,em] = betaNMF(X,r,options);
% Plot 
figure; 
plot(es); hold on; plot(em,'--'); 
h = legend('MU', 'Modified MU'); 
set(h,'Interpreter','latex'); 
xlabel('Iterations','Interpreter','latex'); 
ylabel('Error: $\frac{}{} \|X-WH\|_F$', 'Interpreter', 'latex'); 
axis([ 30 100 em(100) es(30) ]); 