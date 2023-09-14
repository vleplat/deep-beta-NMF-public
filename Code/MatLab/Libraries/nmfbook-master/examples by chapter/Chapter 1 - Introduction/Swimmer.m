% Different algorithms applied on the Swimmer data set
clear all; clc; 

load Swimmer; 
X = X'; 
[m,n] = size(X); 
r = 17; 
%% Standard NMF fails to recover the parts (although the error is 0) 
[Wnmf,Hnmf] = FroNMF(X,r); 
affichage(Hnmf',17,20,11); 
title('Basis images with standard NMF', 'Interpreter', 'latex');
%% Four methods that are able to perfectly decompose the Swimmer data set; 
%% see Chapter 8 - NMF algorithms
%% 1. Separable NMF 
[K,Hsnpa] = SNPA(X,r); 
affichage(Hsnpa',17,20,11); 
title('Basis images with separable NMF', 'Interpreter', 'latex');
% You can obtain the same results as with separable NMF with ONMF, 
% min-vol NMF and NMU, if properly initialized.  
%% 2. ONMF -- very sensitive to initilization / default = SNPA
[Wonmf,Honmf] = alternatingONMF(X,r);
affichage(Honmf',17,20,11); 
title('Basis images with ONMF', 'Interpreter', 'latex');
%% 3. min-vol NMF
[Wminvol,Hminvol] = minvolNMF(X,r);
affichage(Hminvol',17,20,11); 
title('Basis images with min-vol NMF', 'Interpreter', 'latex'); 
%% 4. NMU 
[Wnmu,Hnmu] = recursiveNMU(X,r); 
Hnmu = Hnmu'; 
affichage(Hnmu',17,20,11); 
title('Basis images with NMU', 'Interpreter', 'latex'); 