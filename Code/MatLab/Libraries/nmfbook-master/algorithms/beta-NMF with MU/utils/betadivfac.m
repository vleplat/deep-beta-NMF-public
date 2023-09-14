% Computation of the beta-divergence between X and W*H, that is, 
% dist = D_beta(X,W*H). This can be done more efficiently that with the
% function betadiv if X is sparse for the KL divergence and the Frobenius
% norm. 

function dist = betadivfac(X,W,H,beta); 

if beta == 0 % Itakura–Saito distance
    XsY = X./(W*H+eps); 
    Z = XsY - log(XsY+eps) - 1; 
    dist = sum( Z(:) ); 
    
elseif beta == 1 % Kullback-Leibler divergence
    if issparse(X)
        patX = (X>0); 
        patWH = blockrecursivecompwiseprodsparselowrank(patX,W,H,@xy); 
        dist =  sum( X(patX).* log( X(patX)./(patWH(patX)+eps)+eps ) ) - sum(X(:)) + sum(sum( H.*(repmat(sum(W)',1,size(X,2))) )) ;  
    else
        Y = W*H; 
        Z = X.*log(X./(Y+eps) + eps) - X + Y; 
        dist = sum( Z(:) ); 
    end
    
elseif beta == 2 % Frobenius norm
    dist = 0.5*(norm(X,'fro')^2 - 2*sum(sum( (W'*X).*H ) ) + sum(sum( (W'*W).*(H*H') ) )); 
    
else % Other beta divergences
    Y = W*H; 
    Z = ( X.^beta + (beta-1)*(Y.^beta) - beta * X.*(Y.^(beta-1)) ) /beta/(beta-1); 
    dist = sum( Z(:) ); 
end