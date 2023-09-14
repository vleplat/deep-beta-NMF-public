% Test whether the necessary condition for SSC1 is easily satisfied on
% randomly generated matrices; see Figure 4.6 in the book. 

clear all; clc; close all; 

% Tests on a r-by-n matrix H 
% with r in the set 
r = [10:20:50];              % in the book: 10:10:50 and 10:10:100 
% with density in the set 
densityofH = [0.1:0.2:0.9];  % in the book: 0.1:0.1:0.9
n = 100;                     % in the book: 100 and 200
nattemps = 2;                % in the book: 100

disp('****************************'); 
fprintf('     Test for n = %2.0f   \n', n); 
disp('****************************'); 
nr = length(r);
nd = length(densityofH);
SSCnec = zeros(nr,nd);
for ir = 1 : nr
    for id = 1 : nd
        for i = 1 : nattemps
            % Generate H randomly
            H = full( sprand(r(ir),n,densityofH(id)) );
            % Resample zero columns
            I = find( sum(H) == 0); 
            while ~isempty(I)
                H(:,I) = full( sprand(r(ir),length(I),densityofH(id)) );
                I = I( find( sum(H(:,I)) == 0) ); 
            end
            % Check the necessary condition for SSC1 to be satisfied
            SSCnecirid = SSC1_nec_cond(H);
            SSCnec(ir,id) = SSCnec(ir,id) + SSCnecirid;
        end
        fprintf('rank = %2.0f, density = %1.1f, number of success = %1.0f/%1.0f. \n', ... 
            r(ir), densityofH(id), SSCnec(ir,id), nattemps)
    end
end
SSCnec = SSCnec/nattemps; 
disp('****************************');  

% Display results 
set(0, 'DefaultAxesFontSize', 25);
figure; 
imagesc(SSCnec(end:-1:1,:)); 
hold on; 
colormap(gray); 
colorbar; 
title(sprintf('$n = %2.0f$', n),'Interpreter','latex'); 
set(gca,'ytick',1:length(r)); 
set(gca,'yticklabel',r(end:-1:1));
ylabel('rank ($r$)','Interpreter','latex');
set(gca,'xtick',1:9); 
set(gca,'xticklabel',densityofH);
xlabel('density ($d$)','Interpreter','latex')