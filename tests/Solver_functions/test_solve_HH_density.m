function [pass, name] = test_solve_HH_density()
% test functionality of eval_sol_freqs
% this will NOT catch whether the auto sampler is "wrong" in some way, just
% checks for bugs that break the code and throw an error

%% 2D Tsunami no paneling
% Incident wave

inobj.wave_type = {'Tsunami', 9, 1, .88, 1+.75*1i, 20};
dimension = 2;

inobj.comp_type_freq = {'known', [], []};

%% obtaining s
inobj.bnd_type = {'Kite'} ;
s = setup_boundary(inobj.bnd_type); % may require CHEBFUN ON PATH.

[~, Uinc, ~, ~, ~, ~] = setup_incident(inobj.wave_type, inobj.comp_type_freq, dimension, s);

M = 75;

% Solver type
inobj.solver_type = 'basic';

j = 1;
name{j} = '2D with Tsunami wave using basic solver without paneling';
pass(j) = 1;
try
    solve_HH_density(Uinc{1}, s, [1,2], linspace(1,2,10), M, inobj.solver_type);
catch ME
    pass(j) = 0;
end


%% 2D Tsunami with paneling
% Incident wave

inobj.wave_type = {'Tsunami', 9, 1, .88, 1+.75*1i, 20};
dimension = 2;

inobj.comp_type_freq = {'known', [], []};

%% obtaining s
inobj.bnd_type = {'Kite'} ;
s = setup_boundary(inobj.bnd_type); % may require CHEBFUN ON PATH.

[~, Uinc, ~, ~, ~, ~] = setup_incident(inobj.wave_type, inobj.comp_type_freq, dimension,s);

M = {10,2};

% Solver type
inobj.solver_type = 'basic';

j = j+1;
name{j} = '2D with Tsunami wave using basic solver with paneling';
pass(j) = 1;
try
    solve_HH_density(Uinc{1}, s, [1,2], linspace(1,2,10), M, inobj.solver_type,0);
catch ME
    pass(j) = 0;
end


%% 3D Tsunami
% Incident wave

inobj.wave_type = {'Tsunami', 9, 1, .88, [0, 0, 1.1], 20};
dimension = 3;

inobj.comp_type_freq = {'known', [], []};

%% obtaining s
inobj.bnd_type = {'Torus',.5,1};
s = setup_boundary(inobj.bnd_type); % may require CHEBFUN ON PATH.

[~, Uinc, ~, ~, ~, ~] = setup_incident(inobj.wave_type, inobj.comp_type_freq, dimension, s);

M = 50;

% Solver type
inobj.solver_type = 'basic';

j = j + 1;
name{j} = '3D with Tsunami wave using basic solver';
pass(j) = 1;
try
    solve_HH_density(Uinc{1}, s, [1,2], linspace(1,2,10), M, inobj.solver_type);
catch ME
    pass(j) = 0;
end


