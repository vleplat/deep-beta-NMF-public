% Compute all vertices as the columns of V for the polytope defined by
% { x | e^T x = 1, Ax >= 0 }

function V = computevert(A)

[m,r] = size(A);
T = nchoosek(1:m,r-1);
lT = length(T); 
if lT > 1e6 
    warning('lT is large');
    V = 0; 
    return;
end
V = [];
e = ones(1,r);
b = zeros(r,1); b(r) = 1;
for i = 1 : lT
    % Select the columns s.t. Sx = b
    S = [A(T(i,:),:); e];
    if rank(S) == r
        x = S\b;
        if min(A*x) >= -1e-6 % vertex belongs to the polyhedron
            if size(V,2) == 0
                V = [V x];
            elseif min( sum( (V - repmat(x,1,size(V,2))).^2 ) ) > 1e-6 
                V = [V x]; % avoids adding repeated vertices  
            end
        end
    end
end
