function [ w ] = sinebell(off, stop, pow, tSize)
%INPUTS:
%       off:       Specifies the starting point of the sine-bell in units of pi radians. Common values are 0.0 (for a sine window which starts height at 0.0) and 0.5 (for a cosine window, which starts at height 1.0). The default value is 0.0
%       stop:      Specifies the ending point of the sine-bell in units of pi radians. Common values are 1.0 (for a window which goes to 0.0 height at the last point) and 0.95 (for a window which doesn't go all the way to 0.0). The default value is 1.0
%       pow:       Specifies the exponent of the sine-bell; Non-integer values are allowed. Common values are 1.0 (for ordinary sine-bell) and 2.0 (for squared-bell functions). The default value is 1.0
%       tsize:     Specifies the size of the window

%OUTPUTS:
%       w: sinebell window
x=0:tSize-1;
w=sin(pi*off+pi*(stop-off).*x./(tSize-1)).^pow;
end

