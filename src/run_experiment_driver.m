%EXPERIMENT DRIVER TEST FILE. 
%
%% shows how to run the experiments_driver.m function

clear;
clf

%%
% 1. shows how to load up the input object
% 2. call driver and build the solution object
% 3. some examples on how to plot aspects of the solution object


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. LOAD PREFERENCES INTO AN INPUT OBJECT
% choose the incident wave and scatterer: 
% Boundary domain type: 
%                   2D: {'Kite'};

%                       {'Star',a,w}; a = indent (0<a<1); w = #petals
%                       {'Star',.3,7}

%                       {'Circles',c,R}; c = centers in array; R =
%                       radii in array;
%                       {'Circles',[0 + 0i, 4*(-0.5+1*1i), 5*(1+.5*1i)^(-1)],[2,1,1]}

%                       {'C Curve',a,b,c,d}; a = corner sharpness (bigger =
%                       sharper); b in (0,pi) = 2(pi - channel arch length) (pi = no channel); 
%                       c = fatness (larger = fatter all around);

%                       d = direction (1 = cavity left; exp(1i*(pi-theta)) = rotation theta)  
%                       {'C Curve',3,2.5,1,1} 
%                       {'C Curve',3,2.5,1,exp(-1i*3*pi/4)}

%                       {'Whisper Gallery',a,b,fat}; a = ?; b = ?; fat =
%                       fatness
%                       {'Whisper Gallery',3,6/5*3,1.1}

%                       {'Crescent Gallery',a,b,d}; a = ?; b = ?; d =
%                       distance 
%                       {'Crescent Gallery',.1,.9,4}

%                       {'Teardrop', m}; m = sharpness
%                       {'Teardrop', 1}

%                       {'Keyhole',r,R,e,tt}; r = small radius; R = large radius;
%                                          e = half tube height; tt= angle
%                                          of hole
%                       {'Keyhole',2,3,.3,pi}


%                       {'Radiator',r,R,e,rout,deltatheta,n,theta}  
%                       {'Radiator',1.5,5,.5,.8,.15,7, pi/3}  

%                       {'Spikes',a}; a+1 = amount of spikes

%                       {'Von Koch Snowflake', a}; a = #iterations
%                       {'Von Koch Snowflake', 3}

%                       {'Perforated Wall',h,w,theta,N}; h = heigth, w =
%                       width, theta = angle; N = #holes

%                       {'"Custom".mat'}; Custom .mat with chebfun func: bd
%                       {'CnC.mat'}
%                       {'pumpkin.mat'}

%                   3D: {'Torus',r,R}; r = smaller radius; R = bigger
%                        radius
%                       {'Torus',.5,1}   

%                       {'Wobbly torus',a,b,c}; a = #petals; b = twist number; c = ? 
%                       {'Wobbly torus',3,2,.2} (change domain size)

%                       {'C Torus',a,b,c,d,ee}
%                       {'C Torus',2.5,2.5,.5,exp(1i*pi/4),2}

% FOR 3D: need to add proxy_LU and fmm (fmm3D is installed)
% inobj.bnd_type = {'C Curve',3,3*pi/4,1,exp(-1i*3*pi/4)};
inobj.bnd_type = {'C Torus',2.5,2.5,.5,exp(1i*pi/4),2};    

%% Incident Wave type: 
% Appropriate freq band is auto-selected. 
% Add: final times

% {'Tsunami', base frequency, wave speed, wave width, wave direction, time offset @ z=0}
% {'Tsunami', 9, 1, .88, 1+.75*1i, 20} % 2D
% {'Tsunami', 9, 1, .88, [0, 0, 1], 20} % 3D

% {'Pulse', base frequency, [centers], pulse width, time offset} 
% Add: pulse speed
% {'Pulse', 15, [2-2i], .6, 12} 

% {'Bump', base frequency, Bump speed, Bump width, wave direction, Bump center @ T, T}
% {'Bump', 9, 1, .88, -1-1i, 7.5+7.5i, 20} % 2D
% {'Bump', 9, 1, .88, [-1; -1; -1], [9; 9 ;9], 20} % 3D

% {'Chirp', [head frequency, tail frequency], wave speed, wave width, wave direction, time offset @ z=0}
% {'Chirp', [9,20], 1, .88, 1+.5*1i, 15} % 2D
% {'Chirp', [9,20], 1, .88, [1, .5, 1.5], 15} % 3D

inobj.wave_type = {'Tsunami', 9, 1, .88, 0+1*1i, 20} ;


%% Type of field to solve for (in time): 'total' = total field, 'scatter' = scattered field only
% 'incident' = just the incident wave (this is overkill for the task)
inobj.field_type = 'total';


%% domain sample points: 
% automatically sample domain or input where in space we should evaluate:
n = 50; 
inobj.autosample = {1, n}; % 1 = call autosampler, 0 = do not call
inobj.dom_range = {[-5,5],[-5,5],[-5,5]}; % Depends on dimensionality. 
% inobj.dom_range = {[-5,5],[-5,5]}; % Depends on dimensionality. 
% if called, points will be subsampled from an n^2 equi 2D grid over domain
inobj.zz = []; % if autosampler is off, give us evaluation points. ow zz = []; 
inobj.domplotme = 1; % plot the sampled points in domain? 


%% time locations: you can add a vector of times to evaluate at. 
% If left empty, times will be auto-assigned based on inc_wave
inobj.times = []; 


%% Fourier transform : {'known'} = Theoretical Fourier transform is known
%                   : {'basic',N_ft} = Fourier transform is computed through num. quad (N points)    
%                   : {'window all', N_ft, N_win} = Fourier transform is windowed through num. quad (N_ft points) per N_win windows (multiple HH solvers)
%                   : {'window one', N_ft, N_win} = Fourier transform is windowed through num. quad (N_ft points) per N_win windows (one HH solvers)
inobj.comp_type_freq = {'known', [], []};


%% Which solver for HH eqns: 'basic' = Zetatrap pckg. O(M^3) HH solver
inobj.NystM = {50,2}; % M = {number of points in Nystrom discr per panel, Added pannels per corner}. 
inobj.panelplotme = 1;
inobj.solver_type = 'basic';


%% Evaluation type for frequency space: 'basic' or 'fmm'
inobj.eval_type_freq = 'fmm'; 
% if you want to evaluate only some subset of frequencies, put them here:
% otherwise, the evaluated frequences will be those needed to do the
% inverse Fourier transform to compute solution in time. 
inobj.eval_freqs = []; 



%% Evaluation type for time sol: 'GLquad' or 'fastsinc' or 'complexify'
% GLquad = GL quadrature over freqs: N = # quad pts. 
% fastsinc = {'fastsinc', N}, N = bandlimit for trig. poly approx to HH solution.
% TO DO: implement happycheck and auto-sampling for this purpose. 
% complexify = {'complexify', N, del}

N = 750;
inobj.eval_type_time = {'GLquad',N, .01};

%% Save solution: 1 or 0
inobj.savesol = 0;
inobj.solname = "Test";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% NOW WE CAN CALL THE DRIVER %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. SOLVE THE PROBLEM/CREATE SOLUTION OBJECT
% NOTE I HAVENT TESTED ALL CASES. WILL BE BUGGY FOR A BIT. 

sol = experiments_driver(inobj); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot Results
% plot: 1 or 0
mkim = 1;
% Save: 1 or 0
vidobj.saveim = 0;

vidobj.im_name = 'C-bump';
vidobj.startpoint = 400;

% plotting options: 'custom', 'positive_negative', 'normal', 'logplot' 
vidobj.plottemplate = 'custom';
if strcmp(vidobj.plottemplate, 'custom')
    
    vidobj.soltype = 'time'; %time/freq

    vidobj.frame = 0; %save specific frame: integer/save video: 0

    vidobj.plottype = {'surf', 2}; %surf or mesh

    vidobj.clr = cmocean('dense'); %color

    vidobj.valuetype = 'real'; %abs, real or imag

    % vidobj.rng = {[xlims], [ylims], [zlims]}
    vidobj.rng = {}; % Automatic range  
    % vidobj.rng = {[-2,2], [-2,2]}; % Automatic z range 
    % vidobj.rng = {[-1,1]};  % Automatic x,y range
    % vidobj.rng = {[-2,2], [-2,2], [-1,1]};  % set ranges

end


if mkim % create a GIF:
    result_plot(sol, vidobj)
end

%% 3. EXAMPLES PLOTS FOR 2D CASE

%%
% plot a slice of the solution for a fixed frequency: 
% k = index number associated with the frequency
%
% for the complexified case (contour + fast sinc, it's a little more
% complicated to do freq plots. We will add this functionality in later.)
% if ~strcmpi(inobj.eval_type_time, 'complexify')
%     freqs = sol.freq; 
%     k = 60; 
%     frq = freqs(k); 
%     style = 'imag'; % style = 'abs', 'real', or 'imag'
%     plottype = {'surf', [.625, 58.09]}; % {'mesh' or 'surf', 2 = flat, 3 = vertical view, [a, b] = view vector}
%     rng = []; %input [xlim, ylim] or leave empty to autoselect
%     clr = cmocean('phase'); %colormap (either the matrix or a string for matlab maps)
%     plot_slice_freq(sol, k,style,plottype, rng, clr)
%     hold on
%     title('Solution at frequency ', frq)
%     hold off
% end
% %%
% % plot a slice of the solution for a fixed time: 
% times = sol.times; 
% k = 1; 
% tm = times(k); 
% style = 'real'; % style = 'abs', 'real', or 'imag'
% plottype = {'surf', 2}; % {'mesh' or 'surf', 2 = flat, 3 = vertical view, [a, b] = view vector}
% rng = []; %input [xlim, ylim] or leave empty to autoselect
% clr = cmocean('topo'); %colormap (either the matrix or a string for matlab maps)
% clf
% plot_slice_time(sol, k,style,plottype, rng, clr)
% %%
% % going through time slices
% clf
% rngz = [-1, 1];
% for j = 1:length(times)
%     plot_slice_time(sol, j,style,plottype, rng, clr)
%     clim('manual')
%     clim(rngz)
%     %zlim('manual')
%     %zlim(rngz)
%     shg
% end