% Applying the lower bounds presented in Chapter 3.3 on the slack matrix of
% the regular octagon 
clear all; clc; 

disp('*****************************************************'); 
disp('Illustration of the lower bounds for the nonnegative'); 
disp('  rank on the slack matrix of the regular octagon'); 
disp('*****************************************************'); 

disp('The octagon can be represented as {x | Ax <= b} with') 
a = sqrt(2)/2;
A = [ 1 a 0 -a -1 -a 0 a;
    0 a 1 a 0 -a -1 -a ]'; 
AT = A' 
b = (1+a)*ones(8,1); 
bT = b' 

disp('Its vertices are the columns of the matrix'); 
V = [(1+a) a;
    a (1+a);
    -a (1+a);
    -(1+a) a;
    -(1+a) -a;
    -a -(1+a);
    a -(1+a);
    (1+a) -a]'

fprintf('Click any button to continue.'); 
pause; 
fprintf('\n');

disp('Its slack matrix is given by'); 
S = b(1) - A*V

fprintf('Click any button to continue.'); 
pause; 
fprintf('\n'); 

disp('An Exact NMF of size 6 is given by S = WH with'); 
W =  [1      0      0      1      0    1+2*a
    1+2*a      0      0      0      1   2+2*a
    1      1      0      0      0   1+2*a
    0      2*a      1      0      0   1
    0      1      1+2*a      0      1  0
    1      0      2+2*a      0      1+2*a  0
    0      0      1+2*a      1      1  0
    0      0      1      2*a      0  1]
H = [0      0      0      1      0      0    1    0
    1      0      0      0      0       1     1+2*a    1+2*a
    1      1      0      0      0      0    0      0
    0      1      1+2*a      1+2*a      1      0    0      0
    0      0      1      0      0      0     0      1
    0      0      0      0      1      1     0      0]

fprintf('We have ||S-WH||_F = %2.2d.\n', norm(S - W*H,'fro')); 

fprintf('Click any button to continue.'); 
pause; 
fprintf('\n'); 

[rc,geo,nnucnorm,tausos,hypsep] = lowerbounds_nnr(S,8); 