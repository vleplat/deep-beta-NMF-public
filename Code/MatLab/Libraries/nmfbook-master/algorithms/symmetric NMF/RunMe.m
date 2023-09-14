% Simple test for our coordinate descent methods on a randomly generated
% matrix. 
% 
% See A. Vandaele, N. Gillis, Q. Lei, K. Zhong and I.S. Dhillon, "Efficient 
% and Non-Convex Coordinate Descent for Symmetric Nonnegative Matrix 
% Factorization", IEEE Trans. on Signal Processing 64 (21), pp. 5571-5584, 
% 2016. 

% Randomly generated symmetric nonnegative matrix A to factorize
n=500; A=rand(n,n);
A=A+A';

% Inner rank of the factorization
r=30;

% Options (see loadoptions.m)
options.maxiter   = 100;
options.timelimit = 5;
options.initmatrix='dense01'; % random initialization
options.seed=0; % Use same seed to have the same initialization 

% Comparison of two CD methods: cyclic updates of the variables vs. 
% random shuffling of the columns. 
options.shuffle_columns = 0;
tic; [H1,e1,t1] = symNMF(A,r,options); toc

options.shuffle_columns = 1;
tic; [H2,e2,t2] = symNMF(A,r,options); toc

plot(t1,e1); hold on; plot(t2,e2,'red');
legend('CD-Cyclic-Rand','CD-Shuffle-Rand');
xlabel('Time (s.)'); 
ylabel('Error  -  1/2 * ||A - HH^T||_F^2'); 