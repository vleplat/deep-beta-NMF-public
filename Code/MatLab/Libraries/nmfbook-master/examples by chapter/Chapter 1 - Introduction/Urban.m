% Hierarchical NMF applied on the Urban hyperspectral image 
% The solution shown in the book (Figure 1.6) is obtained with the 
% hierarchical clustering technique from 
% N. Gillis, D. Kuang and H. Park, "Hierarchical Clustering of 
% Hyperspectral Images using Rank-Two Nonnegative Matrix Factorization", 
% IEEE Trans. on Geoscience and Remote Sensing 53 (4), pp. 2066-2078, 2015.
% It is briefly described in Section 8.8.2 Hierarchical/divide-and-conquer 
% approaches. 
clear all; clc; 

load Urban; 
r = 6; 
% You will have to choose the clusters to split; 
% see also Figure 8.6 to see which cluster to split to obtain the solution
% shown on Figure 1.6
disp('---------------------------------------------------------------------------------'); 
disp('You will need to use the following answers to get the same results as in the book:') 
disp('Do you want to visually choose the cluster to be split? y'); 
disp('What is the number of pixels in each row of your hyperspectral images? 307')
disp('Which cluster do you want to split? 2 - 1 - 3 - 5 - 0'); 
disp('---------------------------------------------------------------------------------'); 
disp('Click any button to start'); 
pause 
[~, W] = hierclust2nmf(X,[],1); 
H = NNLS(W,X); 
affichage(H',3,307,307); 
title('From left to right, top to bottom: grass, trees, roof tops 1, dirt, road, roof tops 2', ... 
'Interpreter', 'latex'); 