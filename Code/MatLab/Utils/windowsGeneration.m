function [windows]= windowsGeneration(winsel,WINDOWSIZE)

addpath('.\Utils');

   switch(winsel)
   case 1 
      fprintf(' ->Hamming windows selected \n' );
      windows=hamming(WINDOWSIZE);
   case 2 
      fprintf(' ->Hann windows selected" \n' );
      windows=hann(WINDOWSIZE);
   case 3 
      fprintf(' ->Blackman windows selected \n' );
      windows=blackman(WINDOWSIZE);
   case 4
      fprintf(' ->Sinebell windows selected \n' );
      [windows]=sinebell(0,1,1,WINDOWSIZE);
   otherwise
      fprintf(' ->Invalid choice\n' );
   end