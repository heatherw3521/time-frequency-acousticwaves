function [pass, name] = test_Clen_Curt()
% test functionality of Clen_Curt
% this will catch whether the result is inaccurate and
% checks for bugs that break the code and throw an error

%%
j = 1; 
pass(j) = 1; 
name{j} = 'adaptive Clenshaw Curtis';

f = @(x) exp(-x.^2);
I = [-10,10];
try 
    Res = Clen_Curt(f,I,0);
catch
    pass(j) = 0;
end
%
j = j+1;
pass(j) = 1; 
name{j} = 'accuracy adaptive Clenshaw Curtis';

if ~exist('Res','var') || abs(Res-1.772453850905516)>1e-12
    pass(j) = 0;
end

j = j+1;
pass(j) = 1; 
name{j} = 'Clenshaw Curtis';

try 
   Res = Clen_Curt(f,I,500);
catch
    pass(j) = 0; 
end
%
j = j+1;
pass(j) = 1; 
name{j} = 'accuracy Clenshaw Curtis';

if ~exist('Res','var') || abs(Res-1.772453850905516)>1e-12
    pass(j) = 0; 
end
%%
j = j+1; 
pass(j) = 1; 
name{j} = 'adaptive Clenshaw Curtis at z points';

f = @(z,x) z .* exp(-x.^2);
I = [-10,10];
try 
    Res = Clen_Curt(f,I,0,linspace(0,1,10));
catch
    pass(j) = 0;
end
%
j = j+1;
pass(j) = 1; 
name{j} = 'accuracy adaptive Clenshaw Curtis at z points';

if ~exist('Res','var') || sum(abs(Res-(linspace(0,1,10).')*1.772453850905516))>1e-12
    pass(j) = 0; 
end
%%
j = j+1; 
pass(j) = 1; 
name{j} = 'Clenshaw Curtis at z points';

f = @(z,x) z .* exp(-x.^2);
I = [-10,10];
try 
    Res = Clen_Curt(f,I,500,linspace(0,1,10));
catch
    pass(j) = 0;
end
%
j = j+1;
pass(j) = 1; 
name{j} = 'accuracy Clenshaw Curtis at z points';

if ~exist('Res','var') || sum(abs(Res-(linspace(0,1,10).')*1.772453850905516))>1e-12
    pass(j) = 0; 
end

