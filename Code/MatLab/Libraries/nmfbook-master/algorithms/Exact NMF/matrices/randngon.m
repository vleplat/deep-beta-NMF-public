% Generates the slack matrix of a generic n-gon (hence the matrix
% has dimension n x n). 
% As n increases it takes more and more time as it becomes more likely for
% a point to belong to the convex hull of the previously generated points. 
% See `Heuristics for Exact Nonnegative Matrix Factorization', A. Vandaele,  
% N. Gillis, F. Glineur and D. Tuyttens, for more details. 

function [S,V,F] = randngon(n)
while 1
    P = randn(n,2);
    o = convhull(P,'simplify',true);
    if length(o) == n+1
        break;
    end
end
V=P-ones(n,1)*mean(P);
for i=1:n
    F(i,1:2) = [V(o(i),2)-V(o(i+1),2) V(o(i+1),1)*V(o(i),1)] / (V(o(i),2)*V(o(i+1),1)-V(o(i),1)*V(o(i+1),2));
    v1=V(o(i),1); v2=V(o(i),2);
    w1=V(o(i+1),1); w2=V(o(i+1),2);
    F(i,1:2) = [v2-w2 w1-v1]/(v2*w1-v1*w2);
end
S=max(ones(n)-F*V',0);