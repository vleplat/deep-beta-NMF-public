function saveEstimates(x_Source,options,myFolder)
disp(' ->Saving data in wav format')
K=options.K;
fs=options.fs;

for i=1:K
   fprintf('  ->Source n° %1.0f saving \n', i); 
   %sound(x_Source(:,i),fs);
   s1 ='source_';
   s2 = num2str(i);
   filename = strcat(myFolder,s1,s2,'.wav');
   audiowrite(filename,x_Source(:,i),fs)
end

end%EOF