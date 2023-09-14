% NMF with the KL divergence applied on the spectrogram corresponding to 
% the dataset `Mary has a little lamb'. 
clear all; clc; 
rng(1234); 
% load data matrix
load piano_Mary
r = 4;  
% KL-NMF; see folder 'algorithms/beta-NMF with MU' 
options.beta = 1; 
fprintf('Computing the solution of KL-NMF. \n'); 
[W,H,e] = betaNMF(X,r,options); 
% Display result 
figure; 
p = 1;
timex = (0:1:(size(H,2)-1))/600*5;
freqx = ((1:size(W,1))*50)/1000; 
notes{1} = 'D_4'; notes{2} = 'C_4'; notes{3} = 'E_4'; notes{4} = 'hammer'; 
for i = 1 : r
    subplot(r,2,p); 
    plot(timex,H(i,:)); grid on;  
    ylabel(notes{i}); 
    if p == 1
        title('Activations of the sources','Interpreter','latex');  
    end
    if i == r
        xlabel('Time (s.)','Interpreter','latex'); 
    end
    subplot(r,2,p+1); 
    plot(freqx,W(:,i)); grid on; 
    ylabel(notes{i}); 
    if p == 1
        title('Frequency response of the sources','Interpreter','latex');      
    end
    axis([ freqx(1) freqx(50) 1e-3 1 ])
    if i == r
        xlabel('Frequency (kHz)','Interpreter','latex'); 
    end
    p = p+2; 
end
figure; 
plot(e); 
xlabel('Iterations','Interpreter','latex'); 
ylabel('$D_{1}(X,WH)$','Interpreter','latex'); 
title('Evolution of the objective function','Interpreter','latex'); 