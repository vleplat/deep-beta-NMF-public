% Perform Y = compwiseprodsparselowrank(X,W,H,beta) 
% by decomposing X in blocks 

function Y = blockrecursivecompwiseprodsparselowrank(X,W,H,fun,nnzparam)

if nargin <= 4
    nnzparam = 1e3; % seems to work well
end

Y = sparse(size(X,1), size(X,2));

% Base step
if nnz(X) < nnzparam
    Y = compwiseprodsparselowrank(X,W,H,fun); 
% Recursive step
else
    nblocks = 2; 
    [m,n] = size(X); 
    mi = ceil(m/nblocks); 
    ni = ceil(n/nblocks); 
    for i = 1 : nblocks
        for j = 1 : nblocks
            indi = 1+(i-1)*mi : min(m,i*mi);
            indj = 1+(j-1)*ni : min(n,j*ni);
            Y(indi,indj) = blockrecursivecompwiseprodsparselowrank(X( indi, indj ),W(indi,:),H(:,indj),fun,nnzparam);
        end
    end
end