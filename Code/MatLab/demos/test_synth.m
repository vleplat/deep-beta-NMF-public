% Test on synthetic data 

close all; clear all; clc; 

rng(2023); 

r = [6 3]; 
n = 1000; 
m = 100; 

W2 = rand(m,r(2)); 
W2 = W2./repmat(sum(W2),m,1); 
H1 = generateH(r(1),1000,3); 
omega = 0.2; % Should be smaller than 0.5 for SSC 
H2 = [omega   1      1   omega   0     0   
        1   omega   0      0   omega   1   
        0      0   omega   1      1   omega]; 
    H2 = H2/(1+omega); % to have sum to one 
W1 = W2*H2; 
X = W1*H1; 
divscal = 1; 

Noise = randn(m,n); 
X=X/divscal + 0.01*Noise/norm(Noise,'fro')*norm(X,'fro');

%% Parameters 
options.maxiter = 200;  % max nu. of it. for init. stage
options.outerit = 300;  % max nu. of it. for our Algorithms
options.min_vol = 1;    % 0: Algorithm-1, 1: Algorithm-2
options.epsi = 10^-11;   % Algorithm 1: 10^-5, Algorithm 2: 10^-4 - 10^-3
%%% min-vol parameters
options.delta = ones(1,length(r));
options.alpha_tilde = [1; 1];
%%%% parameters for ADMM procedure
%%%%% Fast configuration - Algorithm 2 faster and accuracies of ADMM-procedure increase 
%%%%%  along iterations of the global scheme
options.rho = 100;           % 10-100: to fine tune
options.thres = 10^-4;      % stopping criterion for ADMM-procedure
options.innerloop = 1;      % inner loop for Step 1 of ADMM-procedure
options.maxIterADMM = 200;  

%% Run deep MF 
[Wl,Hl,el,inWH,output] = deepKL_NMF(X,r,options);

%%% Errors
figure; 
plot(el./repmat(el(1,:)+1e-16,size(el,1),1))
title('Errors divided by errors at iteration 1'); 
legend('Level 1', 'Level 2', ...   %'Level 3', 
        'Weighted'); 
    
fprintf(' Deep NMF W-Error level 1 = %2.2f%%, level 2 = %2.2f%%\n', 100*compareWs( Wl{1}, W1 ),100*compareWs( Wl{2}, W2 ))
fprintf('Multi NMF W-Error level 1 = %2.2f%%, level 2 = %2.2f%%\n', 100*compareWs( inWH.W{1}, W1 ) ,  100*compareWs( inWH.W{2}, W2 ))

% Figure default
set(0, 'DefaultAxesFontSize', 25);
set(0, 'DefaultLineLineWidth', 2);

X = X*divscal; 
rng(8)
Proj = randn(2,100); 
Xr = Proj*X; 
figure; plot(Xr(1,:),Xr(2,:),'bo'), hold on; 
W1r = Proj*W1; 
plot(W1r(1,:),W1r(2,:),'rx'); 
W2r = Proj*W2; 
plot(W2r(1,:),W2r(2,:),'ks'); 


W1dr = Proj*Wl{1}; 
plot(W1dr(1,:),W1dr(2,:),'go'); 

W2dr = Proj*Wl{2}; 
plot(W2dr(1,:),W2dr(2,:),'yo'); 