function [x_Source,Mask] = temporal_Reconstruction(W,H,Phi,options,V)

%INPUTS:
%       x:                     mixed input signal
%       W:                     Matrix dictionnaries for target spectrogram
%       H:                     Matrix activations for target spectrogram
%       Phi:                   Phase Matrix
%       options.hopsize        Hopsize (number of points between two
%       successive windows)
%       options.WINDOWSIZE     windows size
%       w=options.windows      windows
%       options.RM:            Reconstruction of sources by ISTFT by
%       "Synthesis" (RM=1) or by Wiener Filtering (RM=2)

%OUTPUTS:
%       x_Source:              Matrix with sources in temporal domain (NxK)
%       Mask:                  Order-3 tensor (F*T*K) with masking coefficient

% Loading parameters
RM=options.RM;
[K,~]=size(H);
hop=options.hopsize;
Nfft=options.WINDOWSIZE;
w=options.windows;

addpath('./Utils');
%% Source generation (using Inverse Short Time Fourier Transform)
if(RM==1)
    % % Methode A (synthesis)
    disp(' ->Start Sources generation by synthesis')
    
    for i=1:K
    
        % % Build the spectrogram of each source (K=number of sources)
        
        x_Source_Mod=W(:,i)*H(i,:);
    
        % % Source (modulus)-> Source (imaginary): Z=a+ib = rho *
        %exp(1i(1 for distinction with i for loop)*teta)
        
        x_Source_Imag=x_Source_Mod.*exp(1i*Phi);
    
        % % Inverse Short Time Fourier Transform for each column of
        %x_Source_Imag and unwindowing to build the vector source
        
        x_Source(:,i) = iSTFT(x_Source_Imag, Nfft, w, hop);
        x_Source(:,i)=real(x_Source(:,i)); %Consider only the real part, 
        %normally the inverse of STFT should be pure real number but due to
        %numerical approximations, imaginary part can be persent as well
    
    end

end

if(RM==2)
    % % Methode B (Filtering)
    disp(' ->Start Sources generation by filtering')
    % % Computation of the phase of the input signal x

    for i=1:K
    
        % % Create Mask Filter
        
        Mask(:,:,i) =  (W(:,i)*H(i,:)./(W*H)); %Section 2.3.2 Lefevre PHD
        
        % % Apply the Mask
        
        x_Source_Mod=Mask(:,:,i).*V;

        % % Source (modulus)-> Source (imaginary): Z=a+ib = rho *
        %exp(1i(1 for distinction with i for loop)*teta)
        
        x_Source_Imag=x_Source_Mod.*exp(1i*Phi);
    
        % Inverse Short Time Fourier Transform for each column of
        %x_Source_Imag and unwindowing to build the vector source
        x_Source(:,i) = iSTFT(x_Source_Imag, Nfft, w, hop);
        x_Source(:,i)=real(x_Source(:,i)); 
        
        %Consider only the real part, 
        %normally the inverse of STFT should be pure real number but due to
        %numerical approximations, imaginary part can be persent as well
    
    end
 
end
end%EOF
