This Matlab codes allows you to compute symmetric nonnegative matrix factorizations (symNMF) using an exact coordinate descent method, as described in the paper 

A. Vandaele, N. Gillis, Q. Lei, K. Zhong and I.S. Dhillon, "Efficient and Non-Convex Coordinate Descent for Symmetric Nonnegative Matrix Factorization", IEEE Trans. on Signal Processing 64 (21), pp. 5571-5584, 2016.


You need a C++ compiler to be able to compile mex files: you need to run 

mex -largeArrayDims symNMFmex.cpp 

before being able to use our code (in particular RunMe.m for a simple test on a randomly generated matrix).  


On Windows, you can for example install Microsoft Windows SDK 7.1; see http://nl.mathworks.com/support/compilers/R2012a/win64.html 

Please report any bug to arnaud.vandaele@umons.ac.be. Thanks!

Good luck! 
