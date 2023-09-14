% tested data sets
close all; clear all; clc; 
cd('./');
addpath(genpath('./'));

numexp = 15; % number of tested data sets
namedat{1} = 'NG20';
namedat{2} = 'ng3sim';
namedat{3} = 'classic';
namedat{4} = 'ohscal';
namedat{5} = 'k1b';
namedat{6} = 'hitech';
namedat{7} = 'reviews';
namedat{8} = 'sports';
namedat{9} = 'la1';
namedat{10} = 'la12';
namedat{11} = 'la2';
namedat{12} = 'tr11';
namedat{13} = 'tr23';
namedat{14} = 'tr41';
namedat{15} = 'tr45';
for it = 1 : numexp
    disp('***********************************'); 
    fprintf('          '); disp(namedat{it}); 
    disp('***********************************');   
    if it == 1
        load('./Text data/NG20');
    elseif it == 2
        load('./Text data/ng3sim');
    elseif it == 3
        load('./Text data/classic');
    elseif it == 4
        load('./Text data/ohscal');
    elseif it == 5
        load('./Text data/k1b');
    elseif it == 6
        load('./Text data/hitech');
    elseif it == 7
        load('./Text data/reviews');
    elseif it == 8
        load('./Text data/sports');
    elseif it == 9
        load('./Text data/la1');
    elseif it == 10
        load('./Text data/la12');
    elseif it == 11
        load('./Text data/la2');
    elseif it == 12
        load('./Text data/tr11');
    elseif it == 13
        load('./Text data/tr23');
    elseif it == 14
        load('./Text data/tr41');
    elseif it == 15
        load('./Text data/tr45');
    end
    
    %X = log(1+Xk); 
    X = Xk; 
    
    r(2) = max(classid);
    r(1) = r(2)*2;
    
    %% Multilayer and Deep NMF
    parameters;
    
    %%% Call of methods
    options.rngseed = 1;
    [Wl,Hl,el,inWH,output] = deepKL_NMF(X,r,options);
    
    %% ------------------------------------------------------------------------
    %%% Post-processing
    %%%------------------------------------------------------------------------
    %%% sanity check for normalization constraints
%     if options.min_vol
%         norm(sum(Wl{1},1) - ones(size(sum(Wl{1},1))))
%         norm(sum(Wl{2},1) - ones(size(sum(Wl{2},1))))
%     else
%         norm(sum(Hl{1},2) - ones(size(sum(Hl{1},2))))
%         norm(sum(Hl{2},2) - ones(size(sum(Hl{2},2))))
%     end
    
    %%% Errors
    figure; 
    subplot(1,2,1); title(namedat{it}); 
    plot(el./repmat(el(1,:)+1e-16,size(el,1),1))
    %title('Errors divided by errors at iteration 1');
    legend('Level 1', 'Level 2', ...   %'Level 3',
        'Weighted');
    
    % classification accuracy of multilayer vs. deep
    H = inWH.H{2}*inWH.H{1};
    [~,classMLNMF] = max(H);
    classMLNMF = bestMap(classid,classMLNMF);
    accMLNMF(it) = 100*sum(classMLNMF == classid)/length(classid); 
    fprintf('Classfication accuracy of multilayer NMF is %2.0f%%\n', accMLNMF(it));
    
    H = Hl{2}*Hl{1};
    [~,classDNMF] = max(H);
    classDNMF = bestMap(classid,classDNMF);
    accDNMF(it) = 100*sum(classDNMF == classid)/length(classid); 

    fprintf('Classfication accuracy of deep NMF is %2.0f%%\n', accDNMF(it));
    
    
    % classification accuracy of NMF
    %[W,H,e] = betaNMF(X,r(2));
    [W,H,e] = multilayerKLNMF(X,r(2),options); 
    H = H{1}; W = W{1}; 
    
    [~,classNMF] = max(H);
    classNMF = bestMap(classid,classNMF);
    accNMF(it) = 100*sum(classNMF == classid)/length(classid); 
    fprintf('Classfication accuracy of NMF is %2.0f%%\n', 100*sum(classNMF == classid)/length(classid));
    
    %%% min-vol info
    if options.min_vol
        subplot(1,2,2); 
        semilogy(output.e_m); hold on;
        semilogy(output.logdetEvol);
        legend('$f(W,H)$','$\log \det(W_1^TW_1 + \delta I)$','$\log \det(W_2^TW_2 + \delta I)$','Interpreter','latex')
        grid on;
        
%         fprintf(' ->Final betadivergence: %0.2f \n', el(end));
%         fprintf(' ->Final penalty term: %0.2f \n', output.e_m(end));
%         fprintf(' ->Initial ratio : %0.2f \n', output.ratio(1,:));
%         fprintf(' ->Final ratio : %0.2f \n', output.ratio(2,:));
    end
    
end

acc = [classMLNMF; accDNMF; accNMF]