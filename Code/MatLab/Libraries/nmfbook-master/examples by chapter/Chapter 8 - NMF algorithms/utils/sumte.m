% reshape error plot from 0 to options.timemax, tmf is the target time. 
% In other words, given the errors em at time tm, so that plot(tm,em) shows 
% the error over time, this code generates the interpolated errors at times
% tmf so that plot(tmf,emf) generates a similar curve. This allows to more
% easily average several runs of the same algorithm on the same dataset. 
% 
% Input:  em and tm are vectors of the same dimensions, tmf is a target time
%         so that tm and tmf have a similar range of values 
%           (tm(1) and tmf(1) should be close, and tm(end)>=tmf(end)). 
% Output: emf so that plot(tmf,emf) interpolates plot(em,em) 

function emf = sumte(em,tm,tmf);  

timemax = max(tmf); 
if max(tm) < timemax
    error('max(tm) should be larger than max(tmf).'); 
end
if max(tm) == timemax
    tm(end) = timemax+1e-6;
end
emf(1) = em(1); 
cnt = 1; 
for i = 2 : length(tmf) 
    while tm(cnt) <= tmf(i)
        cnt = cnt+1;
    end
    if cnt == 1 % this means that tm(1) >= tmf(2), we do not have errors beore tmf(2)
        emf(i) = em(1); 
    else
        lambda = (tmf(i) - tm(cnt-1))/(tm(cnt)-tm(cnt-1)); 
        emf(i) = lambda*em(cnt) + (1-lambda)*em(cnt-1);
    end
end
