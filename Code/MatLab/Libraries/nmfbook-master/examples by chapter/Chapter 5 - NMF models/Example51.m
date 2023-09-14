% Fro-NMF vs. IS-NMF under/over-approximating X; 
% see Example 5.1 

clear all; clc; 
% Dimensions and rank
m = 100; 
n = 100; 
r = 10;
% Number of randomly generated matrices
numtest = 10; % in the book, = 100
options.display = 0; 
fprintf('Running %2.0f experiments:\n', numtest); 
for i = 1 : numtest
    fprintf('%2.0f...', i); 
    X = full( sprand(m,n,0.5) ) + 1e-6;
    W0 = rand(m,r);
    H0 = rand(r,n);
    options.W = W0;
    options.H = H0;
    
    options.beta = 0;
    [W0,H0,e0] = betaNMF(X,r,options);
    R0 = X-W0*H0;
    res0(i,:) = [norm(max(R0,0),'fro') norm( max(-R0,0 ),'fro')]/norm(R0,'fro'); 
    options.beta = 2;
    [W2,H2] = betaNMF(X,r,options); 
    % Note, it would be more efficient to use [W2,H2] = FroNMF(X,r,options); 
    R2 = X-W2*H2;
    res2(i,:) = [norm( max(R2,0),'fro') norm( max(-R2,0 ),'fro')]/norm(R2,'fro'); 
    if mod(i,10) == 0
        fprintf('\n'); 
    end
end
fprintf('\n'); 
% Display results 
disp('Solutions (W,H) for IS-NMF satisfy on average:')
fprintf('||max(0,X-WH)||_F/||X-WH||_F <= %2.2f %%\n', 100*mean(res0(:,1))); 
fprintf('||max(0,WH-X)||_F/||X-WH||_F >= %2.2f %%\n', 100*mean(res0(:,2))); 
disp('Solutions (W,H) for Fro-NMF satisfy on average:')
fprintf('||max(0,X-WH)||_F/||X-WH||_F >= %2.2f %%\n', 100*mean(res2(:,1))); 
fprintf('||max(0,WH-X)||_F/||X-WH||_F <= %2.2f %%\n', 100*mean(res2(:,2))); 