% tested data sets 

numexp = 15; % number of tested data sets 

timeinit = cputime; 

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
    if it == 1
        load('C:/Users/Nicolas Gillis/Dropbox/Data/text mining/NG20');
    elseif it == 2
        load('C:/Users/Nicolas Gillis/Dropbox/Data/text mining/ng3sim');     
    elseif it == 3
        load('C:/Users/Nicolas Gillis/Dropbox/Data/text mining/classic');
    elseif it == 4
        load('C:/Users/Nicolas Gillis/Dropbox/Data/text mining/ohscal');
    elseif it == 5
        load('C:/Users/Nicolas Gillis/Dropbox/Data/text mining/k1b');    
    elseif it == 6
        load('C:/Users/Nicolas Gillis/Dropbox/Data/text mining/hitech');
    elseif it == 7
        load('C:/Users/Nicolas Gillis/Dropbox/Data/text mining/reviews');    
    elseif it == 8
        load('C:/Users/Nicolas Gillis/Dropbox/Data/text mining/sports');    
    elseif it == 9
        load('C:/Users/Nicolas Gillis/Dropbox/Data/text mining/la1');
    elseif it == 10
        load('C:/Users/Nicolas Gillis/Dropbox/Data/text mining/la12');    
    elseif it == 11
        load('C:/Users/Nicolas Gillis/Dropbox/Data/text mining/la2');
    elseif it == 12
        load('C:/Users/Nicolas Gillis/Dropbox/Data/text mining/tr11');
    elseif it == 13
        load('C:/Users/Nicolas Gillis/Dropbox/Data/text mining/tr23');    
    elseif it == 14
        load('C:/Users/Nicolas Gillis/Dropbox/Data/text mining/tr41');    
    elseif it == 15
        load('C:/Users/Nicolas Gillis/Dropbox/Data/text mining/tr45');    
    end
    
    X = dtm'; 
    rng(2020); 
    r = 20; 
    [W,H] = FroNMF(X,r);
    % Subampled words (rows of X)  
    K = []; 
    maxowrdspertop = 25; 
    for i = 1 : r
        [~,coli] = sort(W(:,i),'descend'); 
        K = unique([K; coli(1:maxowrdspertop)]); 
    end
    Xk = X(K,:); 
    save(namedat{it}, 'Xk', 'classid'); 
end