function [x,fs]= loadAudioFile(options)

file = options.fileselection;
   switch(file)
   case 1 
      fprintf('->Loading Synthetic Bass and Drum sample \n' );
      [x,fs]=audioread('./sound_Examples/Drums_Bass/Drum+Bass.wav');
   case 2 
      fprintf('->Loading Piano sample "Mary had a little lamb" \n' );
      [x,fs]=audioread('./sound_Examples/Piano/Mary.wav');
   case 3 
      fprintf('->Loading Prelude Jean-Sebastian Bach sample \n' );
      [x,fs]=audioread('./sound_Examples/Piano_Bach/Prelude.wav');
   
   case 99
       fprintf('->Loading my file \n' );
       [x,fs]=audioread('./Mysample/myfile.wav');
       [~,n]=size(x);
       if(n>1)
           x=sum(x,2);
       end
  
   otherwise
      fprintf('->Invalid choice\n' );
   end