function [pass, name] = test_solve_HH_to_wave()
% test functionality of eval_sol_freqs
% this will NOT catch whether the auto sampler is "wrong" in some way, just
% checks for bugs that break the code and throw an error


dimensions = {{2,{'Kite'}},{3,{'Torus',.5,1}}};
eval_type_times = {{'GLquad', 10}, {'fastsinc', 10}, {'complexify', 10, .1}};
comp_type_freqs = {{'known', [], []},{'basic', 10, []},{'window all', 10, 10},{'window one', 10, 10}};

field_type = 'total';
eval_type_freq = 'basic';
solver_type = 'basic';
wave_type = {'Tsunami', 9, 1, .88, 1+.75*1i, 20};
times = [];
freqs = [];
M = 10;

j = 1;
for i1 = 1:length(dimensions)
for i2 = 1:length(comp_type_freqs)
for i3 = 1:length(eval_type_times)

dimension = dimensions{i1}{1};
bnd_type = dimensions{i1}{2};
comp_type_freq = comp_type_freqs{i2};
inobj.eval_type_time = eval_type_times{i3};
eval_type_time =  eval_type_times{i3}{1};

if dimension == 2
  zzk = [-5 - 5i; -5 ;-5 + 5i;-5 - 5i;-5 + 0i;-5 + 5i;-5 - 5i;-5 + 0i;-5 + 5i];
  wave_type = {'Tsunami', 9, 1, .88, 1+.75*1i, 20};
elseif dimension == 3
  zzk = [-5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5    -5;-5     0     5    -5     0     5    -5     0     5    -5     0     5    -5     0     5    -5     0     5    -5     0     5    -5     0     5    -5     0     5;-5    -5    -5    -5    -5    -5    -5    -5    -5     0     0     0     0     0     0     0     0     0     5     5     5     5     5     5     5     5     5];
  wave_type = {'Tsunami', 9, 1, .88, [0, 0, 1], 20};
end

N = inobj.eval_type_time{2}; sol.szIFT = N;
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
sol.NystM = M;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% set up incident way freq (Uinc) and time (uinc) function handles,
% give bands over which they are supported:
s = setup_boundary(bnd_type); % may require CHEBFUN ON PATH.

sol.s = s;
sol.gam = s.Z;

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

sol = solve_HH_eqn(sol,zzk,0);

name{j} = [num2str(dimensions{i1}{1}),'D using ',eval_type_times{i3}{1},' inverse Fourier transform and, ',comp_type_freqs{i2}{1},' Fourier transform'];
pass(j) = 1;

try
    solve_HH_to_wave(sol,zzk); 
catch ME
    pass(j) = 0;
end
j = j+1;
end
end
end