% Computation of the beta-divergence between X and Y, that is, 
% dist = D_beta(X,Y) = sum(Z(:)) where Z is the component-wise divergence 
% between the entries of X and Y. 

function [dist,Z] = betadiv(X,Y,beta); 

if beta == 0 % Itakura–Saito distance 
    XsY = X./(Y+eps); 
    Z = XsY - log(XsY+eps) - 1; 
    dist = sum( Z(:) ); 
elseif beta == 1 % Kullback-Leibler divergence
    Z = X.*log(X./(Y+eps) + eps) - X + Y; 
    dist = sum( Z(:) ); 
else % Other beta divergences
    Z = ( (X+eps).^beta + (beta-1)*(Y.^beta) - beta * X.*(Y.^(beta-1)) ) /beta/(beta-1); 
    dist = sum( Z(:) ); 
end