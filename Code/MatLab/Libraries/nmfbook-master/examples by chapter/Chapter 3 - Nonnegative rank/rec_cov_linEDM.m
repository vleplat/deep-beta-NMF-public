% Computing the rectangle covering of the linear EDM X_6 
clear all; clc; 

tic; 
n = 6; 
for i = 1 : n
    for j = 1 : n
        X(i,j) = (i-j)^2;
    end
end
X 
[rcX,rec] = rec_cov_bound(X); 
fprintf('The rectangle covering bound of X is %2.0f.\n',rcX);
disp('The rectangles covering X are:')
for i = 1 : rcX
    rec{i},
end
toc; 