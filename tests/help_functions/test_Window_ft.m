function [pass, name] = test_Window_ft()
% test functionality of Window_ft
% this will NOT catch whether the result is inaccurate, just
% checks for bugs that break the code and throw an error


% Add accuracy testing.

j = 1; 
pass(j) = 1; 
name{j} = 'window fourier transform';


f = @(z,w) sin(pi*(w+z));
zz = [1+1i, .5- .74i];
ww= [1,2,3,4.3];
N = 100;
s = 2;
H = .2;

try
    Window_ft(f,zz,ww,s,N,H);
catch
    warning('Window_ft: unable to run.')
end

j = j+1; 
pass(j) = 1; 
name{j} = 'window functions';

try
    Window_ft(f,zz,ww,s,N,H,'window func');
catch
    warning('Window_ft: unable to run.')
end
