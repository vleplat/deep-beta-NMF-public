% Geometric bound for the nonnegative rank based on the rank of the input
% matrix (rankX) and its restricted nonnegative rank (restrnnrankX). 
% See the paper 
% Gillis, N., & Glineur, F. (2012). On the geometric interpretation of the 
% nonnegative rank. Linear Algebra and its Applications, 437(11), 2685-2712.
% or Section 3.4.6 for mode details. 
% 
% Input : rankX of a matrix X, and its restricted nonnegative rank restrnnrankX
% Output: a lower bound lowerbnd on the nonnegative rank of X 

function lowerbnd = geometric_bound(rankX,restrnnrankX)

lowerbnd = rankX; 

while phifct(rankX,lowerbnd) < restrnnrankX
    lowerbnd = lowerbnd+1;
end

end

function phirrp = phifct(r,rp)

phirrp = rp; % value of phi for rw=r
for rw = r+1 : rp
    if r > 3
        phirrp = max(phirrp, faces(rp,rw-1,rw-r) );
    else
        % This is an improvement in the rank-three case 
        phirrp = max(phirrp, min( faces(rp,rw-1,rw-r) , faces(rp,rw-1,rw-r+1)));
    end
end

end

% This function computes the number of k-faces of a polytope in dimension d
% with n vertices.
% k and d must be scalars here.
% k, n and d must be nonnegative integers.
function sum = faces(n,d,k)
if mod(d,2) == 0
    sum = 0;
    l = d/2;
    for i = 0:l-1
        sum = sum + (comb(d-i,k+1-i) + comb(i,k+1-d+i))*comb(n-d-1+i,i);
    end
    sum = sum + (comb(d-l,k+1-l) + comb(l,k+1-d+l))*comb(n-d-1+l,l)/2;
else
    sum = 0;
    l = floor(d/2);
    for i = 0:l
        sum = sum + (comb(d-i,k+1-i) + comb(i,k+1-d+i))*comb(n-d-1+i,i);
    end
end
end

% This function compute a!/[(a-b)!b!] in a clever way to avoid overhead
function c = comb(a,b)
if a < 0 || b < 0 || a-b < 0
    c = 0;
else
    if b > a-b
        c = 1;
        t = 1;
        for i = (b+1):1:a
            c = c*i/t;
            t = t + 1;
        end
    else
        c = 1;
        t = 1;
        for i = (a-b+1):1:a
            c = c*i/t;
            t = t + 1;
        end
    end
end
end