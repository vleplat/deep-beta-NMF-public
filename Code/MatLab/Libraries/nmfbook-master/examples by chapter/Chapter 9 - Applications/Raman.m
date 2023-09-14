% Application of NMF for Self-modeling curve resolution with Raman spectra; 
% see Section 9.3 
clear all; clc; 

load RamanSMCR; 

[m,r] = size(W); 
figure; 
De = 280; 
plot(De:1:De+m-1, W); hold on; 
[r,n] = size(H); 

axis([De De+m-1 0 max(W(:))]); 
xlabel('Wavenumber (cm$^{-1}$)', 'Interpreter', 'latex'); 
ylabel('Intensity', 'Interpreter', 'latex'); 
title('Columns of $W$ - spectral signatures of the compounds', 'Interpreter', 'latex'); 
h = legend('A', 'B', 'C', 'D', 'E'); 
set(h,'Interpreter','latex'); 

figure; 
plot(0:1/3:50,H'); 
xlabel('Time (s.)', 'Interpreter', 'latex'); 
ylabel('Concentration', 'Interpreter', 'latex'); 
title('Rows of $H$ - activation over time of the compounds', 'Interpreter', 'latex'); 

axis([0 50 0 1]); 
h = legend('A', 'B', 'C', 'D', 'E');
set(h,'Interpreter','latex'); 