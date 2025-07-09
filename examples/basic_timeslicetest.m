clear all
close all
%{'Radiator',r,R,e,rout,deltatheta,n,theta}
inobj.bnd_type = {'Radiator',3.8,8.5,2,1.5,.25,6, 5*pi/4};
s = setup_boundary({'Radiator',3.8,8.5,2,1.5,.25,6, 5*pi/4});
clf
plot_scatterer3(s.Z, 1, 1)
axis equal
shg
%%
inobj.wave_type = {'Tsunami', 9, 1, .88, -exp(1i*5*pi/4), 30};
%inobj.wave_type = {'Tsunami', 20, 1, .88, -exp(1i*5*pi/4), 30};
%inobj.bnd_type = {'Keyhole',2,3,.3,pi};
inobj.field_type = 'total';
n = 350; 
inobj.autosample = {1, n};  
inobj.dom_range = {[-10,10],[-10,10]};
inobj.zz = [];
inobj.domplotme = 0;
inobj.times = linspace(8, 150, 300); 
inobj.comp_type_freq = {'known', [], []};
%inobj.NystM = {30,2}; 
%inobj.NystM = {30,4};
inobj.NystM ={5e3,2};
inobj.panelplotme = 0;
inobj.solver_type = 'basic';
inobj.eval_type_freq = 'fmm'; 
inobj.eval_freqs = []; 
N = 900;
inobj.eval_type_time = {'complexify', N, .025};
inobj.savesol = 0;
inobj.solname = "radiator_test";
%%
sol = experiments_driver(inobj); 
alldone('codeswork.mat')
%%
% check plot in time
 clf
%rngz = [-1, 1];
rngz = [-1, 1];
gam =sol.gam;
% % plot a slice of the solution for a fixed time: 
 times = sol.times; 
 k = 2; %14 %100
 tm = times(k); 
 style = 'real'; % style = 'abs', 'real', or 'imag'
 plottype = {'surf', 2}; % {'mesh' or 'surf', 2 = flat, 3 = vertical view, [a, b] = view vector}
 rng = []; %input [xlim, ylim] or leave empty to autoselect
 %clr = cmocean('thermal'); %colormap (either the matrix or a string for matlab maps)
 clr = 'jet';
 clf
 plot_slice_time(sol, k,style,plottype, rng, clr, [], 0)
 %clim('manual')
  %  clim(rngz)
 plot_scatterer3(gam, 0, 0)
 axis equal
 %colormap(slanCM('inferno')) %for mag plot
 colormap(slanCM('iceburn'))
 grid off
 axis tight
 axis off
 set(gcf, 'Color', 'w')
 set(gca, 'fontsize', 25)
 shg
 %%