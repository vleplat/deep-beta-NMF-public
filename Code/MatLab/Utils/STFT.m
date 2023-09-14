function X = STFT(x, Nfft, w, hop)

if nargin < 2;  Nfft = 256; end
if nargin < 3;  w = Nfft; end
if nargin < 4;  hop = 0; end

% expect x as a row
if size(x,1) > 1
  x = x';
end

% Particular case where the window is not provided
if length(w) == 1
  if w == 0
    % special case: rectangular window
    win = ones(1,Nfft);
  else
    if rem(w, 2) == 0   % force window to be odd-len
      w = w + 1;
    end
    halflen = (w-1)/2;
    halff = Nfft/2;   % midpoint of win
    halfwin = 0.5 * ( 1 + cos( pi * (0:halflen)/halflen));
    win = zeros(1, Nfft);
    acthalflen = min(halff, halflen);
    win((halff+1):(halff+acthalflen)) = halfwin(1:acthalflen);
    win((halff+1):-1:(halff-acthalflen+2)) = halfwin(1:acthalflen);
  end
else
  win = w;
end

Nw = length(win);
% now can set default hop (50%)
if hop == 0
  hop = floor(Nw/2);
end

Q = floor(Nw/hop);
x = [zeros(1,(Q-1)*hop) x zeros(1,(Q-1)*hop)];
l = length(x);
T = ceil((l-Nw)/hop)+1;
x = [x zeros(1,(T-1)*hop+Nw-l)];

X = zeros(Nfft/2+1,T);

for t = 1:T,
    time = 1 + (t-1) * hop : Nw + (t-1) * hop;
    xw = x(time)' .* win;
    aux = fft(xw,Nfft);
    aux = aux(1:Nfft/2+1);
    X(:,t) = aux;
end

end
