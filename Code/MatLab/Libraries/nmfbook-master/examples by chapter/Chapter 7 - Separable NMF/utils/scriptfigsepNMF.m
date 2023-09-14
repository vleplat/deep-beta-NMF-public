% display results 

function a = scriptfigsepNMF(results, timings, delta, nalgo, xp)

colors{1} = 'd-' ;   % SPA
colors{2} = 'x--' ; % VCA
colors{3} = '*-' ; % FastAnchorWords 
colors{4} = 's:' ; % SNPA
colors{5} = '^-.' ; % MVE-SPA
colors{6} = 'o-' ; % SPA-SPA 
colors{7} = 'p:' ; % FGNSR 
colors{8} = 'h:' ; % SD-LP 

figure; 
allleg = ['SPA    '; 
          'VCA    '; 
          'FAW    ', 
          'SNPA   ', 
          'MVE-SPA', 
          'SPA-SPA', 
          'FGNSR  ',
          'SD-LP  ',
          ]; 
      
legendcap = []; 
% Running time
disp('***********************************');
disp('Total running time: ');  
for i = 1 : length(nalgo) 
    if xp <= 2
        plot(delta, results(i,:), colors{i}, 'LineWidth', 2, 'MarkerSize', 15); hold on;
    else
        semilogx(delta, results(i,:), colors{i}, 'LineWidth', 2, 'MarkerSize', 15); hold on;
    end
    fprintf([allleg(nalgo(i),:) ' = %2.2f, \n'], timings(i)); 
end
xlabel('Noise level ($\delta$)', 'Interpreter','latex', 'FontSize',25); 
ylabel('Average accuracy', 'Interpreter', 'latex', 'FontSize',25); 
leg1 = legend(allleg(nalgo,:));
set(leg1,'Interpreter','latex'); 

axis([min(delta) max(delta) min(results(:)) 1])
grid on; 

disp('***********************************');
% Robustness
disp('Robustness 100%: ');  
for i = 1 : length(nalgo)
    res = results(i,:) >= 1-1e-16; 
    if res == 1
        rob(i) = delta(end);
    else
        [a,b] = min(res);  
        if b == 1, rob(i) = 0; else rob(i) = delta(b-1); end
    end
    if xp <= 2
        fprintf([allleg(nalgo(i),:) ' = %2.2f, \n'], rob(i)); 
    else
        fprintf([allleg(nalgo(i),:) ' = %2.2d, \n'], rob(i)); 
    end 
end
disp('***********************************'); 

% Robustness
disp('Robustness 95%: ');  
for i = 1 : length(nalgo)
    res = results(i,:) >= 0.95; 
    if res == 1
        rob(i) = delta(end);
    else
        [a,b] = min(res);  
        if b == 1, rob(i) = 0; else rob(i) = delta(b-1); end
    end
    if xp <= 2
        fprintf([allleg(nalgo(i),:) ' = %2.2f, \n'], rob(i)); 
    else
        fprintf([allleg(nalgo(i),:) ' = %2.2d, \n'], rob(i)); 
    end 
end
disp('***********************************'); 
a = 1; 