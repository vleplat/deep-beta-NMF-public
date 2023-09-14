% Computing the rectangle covering bound for matrix X 
%
% This implementation is rather naive, and can be improved. 
% 
% Input : - a nonnegative matrix X 
%         - options (optional): see the documentation of the function intlinprog
% Output: - rcX is the rectangle covering number of X 
%         - rec{i} is the ith rectangle needed to cover X in an optimal
%               rectangle covering of X, 1 <= i <= rcX. 

function [rcX,rec] = rec_cov_bound(X,options) 

[m,n] = size(X); 
if min(m,n) > 10
    warning('min(m,n) is rather large... you might need to wait a bit'); 
end
if m > n
    X = X'; 
    [m,n] = size(X); 
end  
% Generate all rectangles covering X (brute-force approach)
I = dec2bin(0:2^m - 1) - '0'; 
I = I(2:end,:); 
A = []; % each column of A will correspond to a rectangle 
for i = 1 : size(I,1)
    % find largest rectangle containing the entries in I(i,:) 
    xpos = find(I(i,:) > 0); 
    if length(xpos)>1
        valcol = min( X(xpos,:) ); 
    else
        valcol = X(xpos,:); 
    end
    ypos = find(valcol > 0); 
    S = zeros(m,n); 
    S(xpos,ypos) = 1; 
    A = [A vec(S)];  
end   
% Solve a linear optimization problem over binary variables
numrec = size(A,2); 
f = ones(numrec,1); % minimize number of rectangles 
% Solve min f^T x s.t. Ax >= b and x binary
binX = vec(X); binX(binX > 0) = 1; 
if nargin <= 1
    chosenrec = intlinprog(f,1:numrec,-A,-binX,[],[],zeros(numrec,1),ones(numrec,1)); 
else
    chosenrec = intlinprog(f,1:numrec,-A,-binX,[],[],zeros(numrec,1),ones(numrec,1),[],options); 
end

rcX = sum(chosenrec); 
indrec = find( abs(chosenrec - 1) < 1e-6);  
for i = 1 : rcX
    rec{i} = reshape(A(:,indrec(i)),m,n);
end