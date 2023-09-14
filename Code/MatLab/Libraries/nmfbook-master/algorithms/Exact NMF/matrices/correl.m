% Generates the submatrix of the slack matrix of the Correlation Polytope.
% Rows and columns are indexed by all vectors in a, b \in {0, 1}^n:
% M_n(a,b) = (1-|a'*b|)^2.
% See `Heuristics for Exact Nonnegative Matrix Factorization', A. Vandaele,  
% N. Gillis, F. Glineur and D. Tuyttens, for more details. 

function S=correl(n)
    M=[];
    for i=1:n
        if i==1
            M = [[M ones(1,1)] ; [M zeros(1,1)]];
        else
            M = [[M ones(size(M,1),1)] ; [M zeros(size(M,1),1)]];
        end
    end
    S=zeros(2^n,2^n);
    for i=1:2^n
        a = M(i,:);
        for j=1:2^n
            b = M(j,:);
            S(i,j)=(1-(sum(a.*b)))^2;
        end
    end
    S=rot90(S,2);
end