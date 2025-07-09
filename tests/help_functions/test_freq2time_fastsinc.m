function [pass, name] = test_freq2time_fastsinc()
% test functionality of freq2time_fastsinc
% this will catch whether the result is inaccurate and
% checks for bugs that break the code and throw an error

j = 1;
pass(j) = 1;
name{j} = 'freq2time_fastsinc';

f = @(x,w) sin(pi*(w+x));
f_ift = @(x,t) -(sin(t) * ((1i * sin(pi * x) * t + pi * cos(pi * x)) * sin(t) + (1i * pi * cos(pi * x) - sin(pi * x) * t) * cos(t))) / (pi * (t^2 - pi^2));
x = randn()+randn*1i;
t = randn();
tol = 1e-7;
try 
   Res = freq2time_fastsinc(f,x,t,[0,2], 10, tol);
catch
    pass(j) = 0;
end

j = j+1;
pass(j) = 1;
name{j} = 'accuracy freq2time_fastsinc';

if ~exist('Res','var') || abs(Res-f_ift(x,t))>tol
    pass(j) = 0;
end

end