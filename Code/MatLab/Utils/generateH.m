% Generate H, r by n, with k non-zeros per columns

function H = generateH(r,n,k) 

H = zeros(r,n); 
for i = 1 : n
    rgi = randperm(r); 
    nonzeros = rgi(1:k); 
    H(nonzeros,i) = sample_dirichlet(ones(k,1), 1);
end