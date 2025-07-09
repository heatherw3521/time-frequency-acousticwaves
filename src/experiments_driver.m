function sol = experiments_driver(inobj)
% Experiment driver
% This code will solve the wave equation on domain of choice, and
% produce various images/videos/data/timing related to solution.
%
% input is an object that specifies various experiment choices.
% THESE ARE DESCRIBED IN DETAIL IN THE run_experiment_driver.m FILE.
%
% output is an object with various properties of interest:
%%%%%% INHERIT FROM INPUT OBJ:
% sol.bnd_type = boundary (scatterer type, string descrip, parameters)
% sol.inc = inc wave type (string descrip)
% sol.field_type = field type (total, scatter or inc)
% sol.solver_type = solver used for HH eqns
% sol.comp_type_freq = eval method to get freq space (e.g., 'window on/all')
% sol.eval_type_freq = eval method used in freq space (e.g., 'fmm')
% sol.eval_type_time = eval method use for inverse FT (e.g., 'fastsinc')
% sol.dom = boundary type (string descrip)
% sol.dom_range = range of domain.
% sol.dompts = spatial points for eval in domain(if supplied)
% sol.times = times we wish to evaluate at
% sol.szIFT = # pts used to discretize IFT itegral
% sol.NystM = disc. size for Nystrom systems.
%

%%%%%%%%%%%%% COMPUTED HERE
% sol.domaingrid = {xx,yy,zz, keep}; (grid info if domain autosampled)
% sol.freq = frequency points where solution has a value.
% sol.tau = density layer functions at each point in sol.freq
% sol.bandlimit = freql (lower and upper limit)
% sol.incf = incident field in time. space
% sol.incF = incident field in freq. space
% sol.gam = boundary as a function handle
% sol.s = quadrature object from zeta-trap

% sol.freqvals = Matrix of sol in freq dom at freqs
% sol.timevals = Matrix of sol in time dom at times
% sol.timsolv = timing for solve HH eqn
% sol.timevalfreq = timing for eval in freq space (not counting solve)
% sol.timeevaltime = timing for eval in time domain
%
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET UP EXPERIMENT:
bnd_type = inobj.bnd_type;
wave_type = inobj.wave_type;
comp_type_freq = inobj.comp_type_freq;
field_type = inobj.field_type;
autosample = inobj.autosample{1};
n = inobj.autosample{2};
zz = inobj.zz;
dimension = 2;
if (bnd_type{1} == "Torus" || bnd_type{1} == "Wobbly torus" || bnd_type{1} == "C Torus")
    dimension = 3;
end
plotme = inobj.domplotme;
times = inobj.times;
solver_type = inobj.solver_type;
if length(inobj.NystM) == 1
    NystM = inobj.NystM{1};
else
    NystM = inobj.NystM;
end
eval_type_freq = inobj.eval_type_freq;
freqs = inobj.eval_freqs;

if isempty(inobj.eval_type_time)
    eval_type_time = 'freq only';
    N =[];
else
    eval_type_time = inobj.eval_type_time{1};
    N = inobj.eval_type_time{2}; sol.szIFT = N;
end
if length(inobj.eval_type_time)==3
    delt = inobj.eval_type_time{3};
else
    delt = [];
end

%Inherit from inobj:
sol.eval_freqs = freqs;
sol.bnd_type = bnd_type;
sol.inc = wave_type;
sol.field_type = field_type;
sol.solver_type = solver_type;
sol.comp_type_freq = comp_type_freq;
sol.eval_type_freq = eval_type_freq;
sol.eval_type_time = eval_type_time;
sol.dom = bnd_type;
sol.szIFT = N;
sol.delt = delt;
sol.NystM = NystM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
s = setup_boundary(bnd_type); % may require CHEBFUN ON PATH.

sol.s = s;
sol.gam = s.Z;

if autosample
    rnge = inobj.dom_range; sol.dom_range = rnge;
    [zz, xx, yy, keep] = autosample_domain(n, rnge, bnd_type, plotme);
    if dimension == 2
        zzk = zz(keep);
    elseif dimension == 3
        zzk = [xx(keep);yy(keep);zz(keep)];
    end
    %return keep mask and orig pts in case they are needed/desired.
else
    zzk = zz; 
    % assumes user didn't give us zz outside domain! (keep not needed).
end

if autosample
    sol.domaingrid = {xx,yy,zz, keep};
else
    sol.domainpts = zz;
end
%auto time based on wave:
%total time = [0, 2*T]
%%
% set up incident way freq (Uinc) and time (uinc) function handles,
% give bands over which they are supported:
[uinc, Uinc, Uinc_slow, T, W, windows] = setup_incident(wave_type, comp_type_freq, dimension, s);

sol.incF = Uinc;
sol.incf = uinc;
sol.incF_slow = Uinc_slow;
sol.windows = windows;
sol.bandlimit = W;
%%
if isempty(times)
    fs = 1/(T);
    tt = (0:fs:T).';
else
    tt = times;
end

sol.times = tt(:);
% tester code to check freq band
%r = 1+.5*1i;
%r = r/abs(r);
%[uinc, Uinc] = setup_tsunami(sol,T, 6, 1/.88,1,r,W, 2, zzk)
%% Solve Helmholtz equation
sol = solve_HH_eqn(sol,zzk, inobj.panelplotme);
%% Solve wave equation
sol = solve_HH_to_wave(sol,zzk);
sol.timevals = {sol.timevals}; % this is just to have consistency with freqvals and how we store solutions.
% if strcmpi(eval_type_time, "complexify")
%     sol.cfreq = sol.freq; 
%     frq = real(sol.freq{2});
%     sol.freq = frq;
% end


if inobj.savesol
    outputDir = fullfile(pwd, 'BroadbandHH Results/Solutions');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir,inobj.solname),'sol')
end

end