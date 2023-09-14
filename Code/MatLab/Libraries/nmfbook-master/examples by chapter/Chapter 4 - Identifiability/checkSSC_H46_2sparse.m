% Example 4.29 - sufficiently scattered condition
clear all; clc;
H = [ 1 1 1 0 0 0
      1 0 0 1 1 0
      0 1 0 1 0 1
      0 0 1 0 1 1]
[ssc1,ssc2,xs,ssc1nec] = isSSC(H);
if ssc1nec == 0
    disp('H does not satisfy the necessary condition for the SSC.'); 
else
    disp('An optimal solution xs of max_x ||x||^2 such that H^T x >= 0, e^T x = 1:')
    xs
    if ssc1 == 1
        fprintf('H satisfies SSC1 because ||xs||^2 = %2.2f <= 1.\n', norm(xs)^2);
    else
        fprintf('H does not satisfy SSC1 because ||xs||^2 = %2.2.f > 1.\n', norm(xs)^2);
    end
    if ssc2 == 1
        disp('H satisfies SSC2 because xs can only be a unit vector.');
    else
        disp('H does not satisfy SSC2 because xs is not a unit vector.');
    end
end