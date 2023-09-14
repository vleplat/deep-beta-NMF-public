% Perform Y = compwiseprodsparselowrank(X,W,H,beta) 
% by decomposing X into blocks 

function Y = blockcompwiseprodsparselowrank(X,W,H,beta,nblocks)

[m,n] = size(X); 
mi = ceil(m/nblocks); 
ni = ceil(n/nblocks); 
Y = sparse(size(X,1), size(X,2));
for i = 1 : nblocks
    for j = 1 : nblocks
        indi = 1+(i-1)*mi : min(m,i*mi); 
        indj = 1+(j-1)*ni : min(n,j*ni); 
        Y(indi,indj) = compwiseprodsparselowrank(X( indi, indj ),W(indi,:),H(:,indj),beta); 
    end
end