% Display 
set(0, 'DefaultAxesFontSize', 25);
set(0, 'DefaultLineLineWidth', 2);
figure; 

numarker = 20; 
ebesta = Inf; 

len = length(tm); 
step = max(1,floor(len/numarker)); 
steps = [1:step:len len]; 
semilogy(tm,em-ebest,'d-','MarkerSize', 15); hold on; 
ebesta = min(ebesta, min(em-ebest)); 

len = length(tma); 
step = max(1,floor(len/numarker)); 
steps = [1:step:len len]; 
semilogy(tma,ema-ebest,'d--','MarkerSize', 15); hold on; 
ebesta = min(ebesta, min(ema-ebest)); 

len = length(ta); 
step = max(1,floor(len/numarker)); 
steps = [1:step:len len]; 
semilogy(ta,ea-ebest,'x--','MarkerSize', 15); hold on; 
ebesta = min(ebesta, min(ea-ebest)); 

len = length(tls); 
step = max(1,floor(len/numarker)); 
steps = [1:step:len len]; 
semilogy(tls,els-ebest,'*-','MarkerSize', 15); hold on; 
ebesta = min(ebesta, min(els-ebest)); 

len = length(th); 
step = max(1,floor(len/numarker)); 
steps = [1:step:len len]; 
semilogy(th,eh-ebest,'s:','MarkerSize', 15); hold on; 
ebesta = min(ebesta, min(eh-ebest)); 

len = length(the); 
step = max(1,floor(len/numarker)); 
steps = [1:step:len len]; 
semilogy(the,ehe-ebest,'s-','MarkerSize', 15); hold on; 
ebesta = min(ebesta, min(ehe-ebest)); 

len = length(tf); 
step = max(1,floor(len/numarker)); 
steps = [1:step:len len]; 
semilogy(tf,ef-ebest,'^-.','MarkerSize', 15); hold on; 
ebesta = min(ebesta, min(ef-ebest)); 

len = length(tad); 
step = max(1,floor(len/numarker)); 
steps = [1:step:len len]; 
semilogy(tad,ead-ebest,'o-.','MarkerSize', 15); hold on; 
ebesta = min(ebesta, min(ead-ebest)); 

grid on; 

h = legend('MU', 'A-MU', 'ANLS', 'ALS', 'A-HALS', 'E-A-HALS', 'FPGM', 'AO-ADMM');
set(h,'Interpreter','latex'); 

ylabel('$\frac{\|X-WH\|_F}{\|X\|_F} - e_{best}$','Interpreter','latex'); 
xlabel('Time (s.)','Interpreter','latex'); 

axis([0 options.timemax ebesta 0.5])