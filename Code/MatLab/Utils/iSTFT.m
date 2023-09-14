function x = iSTFT(X, Nfft, w, hop)

if nargin < 2; Nfft = 2*(size(X,1)-1); end
if nargin < 3; w = 0; end
if nargin < 4; hop = 0; end  % will become winlen/2 later

s = size(X);
if s(1) ~= (Nfft/2)+1
  error('number of rows should be fftsize/2+1')
end
T = s(2);
 
if length(w) == 1
  if w == 0
    % special case: rectangular window
    win = ones(1,Nfft);
  else
    if rem(w, 2) == 0 
      w = w + 1;
    end
    halflen = (w-1)/2;
    halff = Nfft/2;
    halfwin = 0.5 * ( 1 + cos( pi * (0:halflen)/halflen));
    win = zeros(1, Nfft);
    acthalflen = min(halff, halflen);
    win((halff+1):(halff+acthalflen)) = halfwin(1:acthalflen);
    win((halff+1):-1:(halff-acthalflen+2)) = halfwin(1:acthalflen);
    win = 2/3*win;
  end
else
  win = w;
end

Nw = length(win);
% now can set default hop
if hop == 0 
  hop = floor(Nw/2);
end

Q=floor(Nw/hop);
xlength = Nw + (T-1)*hop;
x = zeros(xlength,1);

zpf = floor(Nfft/Nw);
Y = X(1:zpf:end,:);

for t = 1:T,
    time = (1 + (t-1) * hop) : (Nw + (t-1) * hop);
    s=ifft(reshape(Y(:,t),[],1),Nw,'symmetric');
    x(time) = x(time) + s;
end;

x = x(((Q-1)*hop+1):(xlength-(Q-1)*hop));

end