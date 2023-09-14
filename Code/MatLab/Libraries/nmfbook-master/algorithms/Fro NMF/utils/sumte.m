% reshape error plot from 0 to options.timemax, 
% tmf is the target time 

function emf = sumte(em,tm,tmf);  

timemax = max(tmf); 
if max(tm) < timemax
    error('max(tm) >= timemax');
end

if max(tm) == timemax
    tm(end) = timemax+1e-6;
end

%tmf = 0:timemax/numarker:timemax; 

emf(1) = em(1); 
cnt = 1; 
for i = 2 : length(tmf) 
    while tm(cnt) <= tmf(i)
        cnt = cnt+1;
    end
    lambda = (tmf(i) - tm(cnt-1))/ (tm(cnt)-tm(cnt-1)); 
    emf(i) = lambda*em(cnt) + (1-lambda)*em(cnt-1);
end