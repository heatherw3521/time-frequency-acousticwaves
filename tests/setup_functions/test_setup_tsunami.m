function [pass, name] = test_setup_tsunami()
% test functionality of plot_setup_tsunami
% this will NOT catch whether the auto sampler is "wrong" in some way, just
% checks for bugs that break the code and throw an err

T = 16;
w0 = 6;
sigma = 1/.88;
c = 1;
W_est = [1,12];

%% 2D
j = 1; pass(j) = 1; name{j} = '2D tsunami';

r = 1-.3i;
bnd_type = {'Kite'};
[zz, xx, yy, ~] = autosample_domain(10, {[-5,5],[-4,5]}, bnd_type, 0);
sol.domaingrid = {xx,yy,zz};
try
    [uinc, Uinc] = setup_tsunami(sol,T, w0, sigma,c,r, W_est, 2, []);
    uinc(2+.3i,1);
    Uinc(2+.3i,1);
catch
    pass(j) = 0;
end

%% 3D
j = j+1; pass(j) = 1; name{j} = '3D tsunami';

r = [1,2,-1];
bnd_type = {'Torus',.5,1};
[zz, xx, yy, ~] = autosample_domain(10, {[-5,5],[-4,5],[-5,5]}, bnd_type, 0);
sol.domaingrid = {xx,yy,zz};
try
    [uinc, Uinc] = setup_tsunami(sol,T, w0, sigma,c,r, W_est, 3, []);
    uinc([1; 2; 3],1);
    Uinc([1; 2; 3],1);
catch
    pass(j) = 0;
end