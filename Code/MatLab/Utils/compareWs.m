function [errW,assignment] = compareWs( W, West )
% Compute     errW = min_perm ||W - West(:,perm)||_F / ||W||_F  (1)
%
% ****** Input ******
% W       : the m-by-r matrix of true basis vectors.
% West    : the m-by-r matrix of estimated basis vectors.
%
% ****** Output ******
% errW    : the relative error given by (1).
%

r = size(W,2); 

% Compute the Euclidean distances between each estimated basis vector and
% each true basis vector
for i = 1 : r
    W(:,i) = W(:,i)/sum(W(:,i));
    West(:,i) = West(:,i)/sum(West(:,i));
end
for i = 1 : r
    for j = 1 : r
        %alpha(i,j) = West(:,j)'*W(:,i) / norm(West(:,j))^2; 
        Dist(i,j) = norm( W(:,i) - West(:,j) , 'fro' )^2; 
    end
end
 
% Assign each estimated basis vector to the closest true one, according to
% the Hungarian algorithm
[assignment,cost] = munkres(Dist); 
% Compute the relative error (1)
errW = norm(W-West(:,assignment),'fro') / norm(W,'fro'); 