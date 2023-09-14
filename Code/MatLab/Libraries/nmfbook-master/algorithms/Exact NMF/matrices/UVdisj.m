% Generates the unique disjointness matrix of order n (hence the matrix
% has dimension 2^n x 2^n). 
% See `Heuristics for Exact Nonnegative Matrix Factorization', A. Vandaele,  
% N. Gillis, F. Glineur and D. Tuyttens, for more details. 

function [S,r]=UVdisj(n)
    if n<=10
        S=1;
        for i=1:n
            if mod(i,2)==0
                S   = [3*S_2 S_2 S_2 S_2 ; S_2 zeros(size(S_2)) S_2 zeros(size(S_2)) ; S_2 S_2 zeros(size(S_2)) zeros(size(S_2)) ; S_2 zeros(size(S_2)) zeros(size(S_2)) S_2];
            else
                S_2 = S;
                S   = [S S ; S zeros(size(S))];
            end
        end
        if nargout>1
            r=nnrank(n);
        end
    else
        display('n must be <=10 (or uncomment the if condition)');
    end
end

function r=nnrank(n)
    if mod(n,2)==0
        r = sqrt(3)^n;
    else
        r = 2*(sqrt(3)^(n-1));
    end
    r = round(r);
end