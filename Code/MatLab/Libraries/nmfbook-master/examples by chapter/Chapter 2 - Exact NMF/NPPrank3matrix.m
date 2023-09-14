% 2-D Representation of a NPP instance corresponding to the RE-NMF 
% instance of a rank-3 nonnegative matrix X; see Theorem 2.11. 
%
% [P,U,V] = NPPrank3matrix(X)
%
% Input.
%   X              : (m x n) matrix or rank 3.
% Output.
%   P (2 x p<=m)         :  vertices of the outer polygon B={x|Fx+g>=0}.
%   U (m x 3), V (3 x n) :  X = UV and columns of U and V sum to one. 
% 
% Code adapted from N. Gillis and F. Glineur, "On the Geometric 
% Interpretation of the Nonnegative Rank", Linear Algebra and its 
% Applications 437 (11), pp. 2685-2712, 2012. 

function [P,U,V] = NPPrank3matrix(X)

if rank(X) ~= 3 || min(X(:)) < 0
    disp('The matrix is not rank 3 or not nonnegative'); 
    return;
end
r = 3; 
[m,n] = size(X);
% 1. Remove zero row/columns and normalize the columns of M
X = X(sum(X')>0,sum(X)>0); 
D = diag(1./sum(X)); X = X*D;
% 2. Compute the decompositon X = UV, U = X(:,K) 
options.display = 0; 
K = SPA(X,3,options); % This selects 3 lin. ind. columns of X, SPA is like 
                      % a modified QR algorithm; see Chapter 7. 
U = X(:,K); 
V = U\X; 
% 3. Find the inequalities of the set P = { x in R^{k-1} | Fx+g >= 0 }
F = U(:,1:r-1)-repmat(U(:,r),1,r-1); 
g = U(:,r);
% 4. Draw P (outer polytope) and S (inner polytope)
P = vertices(F,g); 
Kp = convhull(P(1,:),P(2,:));
figure; plot(P(1,Kp),P(2,Kp),'rx-','MarkerSize',20); 
hold on;
Kv = convhull(V(1,:),V(2,:)); 
plot(V(1,Kv),V(2,Kv),'bo-.','MarkerSize',15);
legend('Outer polygon $\mathcal{B}$', 'Inner polygon $\mathcal{A}$', ... 
    'Interpreter','latex'); 
grid on; 
