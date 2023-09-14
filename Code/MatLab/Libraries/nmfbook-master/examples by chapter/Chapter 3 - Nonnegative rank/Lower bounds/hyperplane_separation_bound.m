% Given Z (m-by-n matrix), this function computes the value of 
% alpha = max_{R in {0,1}^{m,n}} <Z,R> s.t. rank(R) = 1.  
% 
% Because of the combinatorial nature of the problem, we do not recommend
% using this code for min(m,n) > 14. 
%
% Input : matrix Z 
% Output: alpha, 
%         and binary solution (xopt,yopt) with R(i,j) = xopt(i)*yopt(j) 
%         for all i,j corresponding to the optimal solution of the problem 
%         above. 
% 
% See the paper 
% Rothvoß, T. (2017). The matching polytope has exponential extension 
% complexity. Journal of the ACM (JACM), 64(6), 41.
% or Chapter 3.3.7 of the book for more details. 

function [alphaZ,xopt,yopt] = hyperplane_separation_bound(Z) 

[m,n] = size(Z); 
if min(m,n) > 14
    warning('min(m,n) is rather large... you might need to wait a bit'); 
end
if m > n
    Z = Z'; 
    [m,n] = size(Z); 
end
alphaZ = -Inf;   
% Try all possible binary rectangles
X = dec2bin(0:2^m - 1) - '0'; 
Y = zeros(2^m, n); 
for i = 1 : 2^m
    xpos = find(X(i,:) > 0); 
    if length(xpos)>1
        valcol = sum( Z(xpos,:) ); 
    else
        valcol = Z(xpos,:); 
    end
    ypos = find(valcol > 0); 
    Y(i,:) = zeros(1,n); 
    Y(i,ypos) = 1; 
    vali = sum( valcol(ypos) ); 
    if vali > alphaZ
        alphaZ = vali; 
        iopt = i;
        xopt = X(i,:); 
        yopt = Y(i,:); 
    end
end 