% Small numerical test for septrisymNMF
% 
% Given A = W S W^T where W is separable, it recovers W and S. 
clear all; clc; 

% Geneate a synthetic W and S, W separable 
r = 5; 
n = 20; 
W = [eye(r); rand(n,r)]; 
for i = 1 : r
    W(:,i) = W(:,i)/sum(W(:,i));
end
disp('True W and S so that A = WSW^T:')
W = W(randperm(n+r),:)
S = rand(r); 
S = (S+S')/2 
A = W*S*W'; 
% Compute a tri-symNMF
disp('Recovered Wt and St so that A = Wt St Wt^T:')
[Wt,St] = septrisymNMF(A,r) % Code in Chapter 8 - NMF algorithms\separable NMF\separable trisymNMF
fprintf('The recovered (Wt,St) satisfy ||A-Wt*St*Wt^T||_F = %2.2d.\n', norm(A - Wt*St*Wt', 'fro')); 