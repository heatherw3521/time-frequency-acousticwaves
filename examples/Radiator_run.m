
%% MADE WITH ADAPTIVE RATIO PANELS

% NOT THE ONE WE USED FOR EXPERIMENTS: 
% see basic_timeslice.m
%%%%%
%%
clear;
clf
inobj.bnd_type = {'Radiator',3.8,8.5,2,1.5,.25,6, 5*pi/4};
s = setup_boundary({'Radiator',3.8,8.5,2,1.5,.25,6, 5*pi/4});
%inobj.bnd_type = {'Radiator',1.5,5,.8,.8,.15,6, 5*pi/4}; %orig version
%inobj.bnd_type = {'Radiator',3.8,8.5,2,1.5,.25,7, 5*pi/4}; %ofat vers
%s = setup_boundary(inobj.bnd_type);
clf
plot_scatterer3(s.Z, 0, 1)

axis equal
shg


%%
inobj.wave_type = {'Tsunami', 9, 1, .88, -exp(1i*5*pi/4), 20};
inobj.field_type = 'total';
n = 350; 
inobj.autosample = {1, n}; 
inobj.dom_range = {[-8,8],[-8,8]};
inobj.zz = []; 
inobj.domplotme = 0; % plot the sampled points in domain? 
inobj.times = linspace(8, 100, 140); 
inobj.comp_type_freq = {'known', [], []};
inobj.NystM = {5e3,2}; % M = {number of points in Nystrom discr total, Added pannels per corner}. 
inobj.panelplotme = 0;
inobj.solver_type = 'basic';
inobj.eval_type_freq = 'fmm'; 
inobj.eval_freqs = []; 
N = 750;
inobj.eval_type_time = {'complexify',N, .01};
%inobj.eval_type_time = {'GLquad',N, .01};
inobj.savesol = 0;
inobj.solname = "Test";
%%
sol = experiments_driver(inobj); 

%% Plot Results
% this makes a GIF
% % plot: 1 or 0
% mkim = 1;
% % Save: 1 or 0
% vidobj.saveim = 1;
% 
% vidobj.im_name = 'total';
% vidobj.startpoint = 400;
% 
% % plotting options: 'custom', 'positive_negative', 'normal', 'logplot' 
% vidobj.plottemplate = 'custom';
% if strcmp(vidobj.plottemplate, 'custom')
% 
%     vidobj.soltype = 'time'; %time/freq
% 
%     vidobj.frame = 0; %save specific frame: integer/save video: 0
% 
%     vidobj.plottype = {'surf', 2}; %surf or mesh
% 
%     vidobj.clr = cmocean('dense'); %color
% 
%     vidobj.valuetype = 'real'; %abs, real or imag
% 
%     % vidobj.rng = {[xlims], [ylims], [zlims]}
%     vidobj.rng = {}; % Automatic range  
%     % vidobj.rng = {[-2,2], [-2,2]}; % Automatic z range 
%     % vidobj.rng = {[-1,1]};  % Automatic x,y range
%     % vidobj.rng = {[-2,2], [-2,2], [-1,1]};  % set ranges
% 
% end
% 
% 
% if mkim % create a GIF:
%     result_plot(sol, vidobj)
% end
% plot_scatterer3(sol.gam, .6, 1);
% title('');xlabel('');ylabel('');colorbar('off'); set(gca,'xtick',[]); set(gca,'xticklabel',[]); set(gca,'ytick',[]); set(gca,'yticklabel',[])

%%
% check plot in time
clf
%rngz = [-1, 1];
rngz = [-1, 1];
gam =sol.gam;
% % plot a slice of the solution for a fixed time: 
 times = sol.times; 
 k = 37; %14 %100

 tm = times(k); 
 style = 'abs'; % style = 'abs', 'real', or 'imag'
 plottype = {'surf', 2}; % {'mesh' or 'surf', 2 = flat, 3 = vertical view, [a, b] = view vector}
 rng = []; %input [xlim, ylim] or leave empty to autoselect
 %clr = cmocean('thermal'); %colormap (either the matrix or a string for matlab maps)
 clr = 'jet';
 clf
 plot_slice_time(sol, k,style,plottype, rng, clr, [], 0)
 %clim('manual')
  %  clim(rngz)
 plot_scatterer3(gam, .82, 0)
 axis equal
 colormap(slanCM('inferno')) %for mag plot
 %colormap(slanCM('iceburn'))
 grid off
 axis tight
 axis off
 set(gcf, 'Color', 'w')
 set(gca, 'fontsize', 25)
 shg