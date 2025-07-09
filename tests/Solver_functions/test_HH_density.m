function [pass, name] = test_HH_density()
% test functionality of eval_sol_freqs
% this will NOT catch whether the auto sampler is "wrong" in some way, just
% checks for bugs that break the code and throw an error.

wave_types = {{'Tsunami', 9, 1, .88, 1+.75*1i, 20}};
dimensions = {{2,{'Kite'}}};
comp_type_freqs = {{'known', [], []},{'basic', 10, []},{'window one', 10, 5},{'window all', 10, 5}};
solver_types = {'basic'};

j = 1;

for i1 = 1:length(dimensions)
for i0 = 1:length(solver_types)
for i2 = 1:length(comp_type_freqs)
for i3 = 1:length(wave_types)

inobj.wave_type = wave_types{i3};
dimension = dimensions{i1}{1};
inobj.bnd_type = dimensions{i1}{2} ;
inobj.comp_type_freq = comp_type_freqs{i2};
inobj.solver_type = solver_types{i0};

s = setup_boundary(inobj.bnd_type); % may require CHEBFUN ON PATH.

inobj.s = s;
inobj.gam = s.Z;

[~, Uinc, ~, ~, ~, ~] = setup_incident(inobj.wave_type, inobj.comp_type_freq, dimension, s);

B = @(z,w) -Uinc{1}(z,w);

inobj.NystM = 75;
M = 10;

s = setupquad(s, M);

inobj.s = s;

name{j} = [num2str(dimension),'D with ', wave_types{i3}{1} ,' wave using ', solver_types{i0} ,' solver and ', comp_type_freqs{i2}{1} ,' Fourier transform without paneling'];
pass(j) = 1;
try
    HH_density(2, B,s, inobj.solver_type);
catch ME
    pass(j) = 0;
end
j = j+1;
end
end
end
end

i1 = 1;
for i0 = 1:length(solver_types)
for i2 = 1:length(comp_type_freqs)
for i3 = 1:length(wave_types)

inobj.wave_type = wave_types{i3};
dimension = dimensions{i1}{1};
inobj.bnd_type = dimensions{i1}{2} ;
inobj.comp_type_freq = comp_type_freqs{i2};
inobj.solver_type = solver_types{i0};
inobj.plot_panel = 0;


[~, Uinc, ~, ~, ~, ~] = setup_incident(inobj.wave_type, inobj.comp_type_freq, dimension,s);

B = @(z,w) -Uinc{1}(z,w);

s = setup_boundary(inobj.bnd_type); % may require CHEBFUN ON PATH.

inobj.s = s;
inobj.gam = s.Z;

inobj.NystM = 75;
M = {10,2};

s = ns_setupquad(s, M, 0);

inobj.s = s;

name{j} = [num2str(dimension),'D with ', wave_types{i3}{1} ,' wave using ', solver_types{i0} ,' solver and ', comp_type_freqs{i2}{1} ,' Fourier transform with paneling'];
pass(j) = 1;
try
    HH_density(2, B,s, inobj.solver_type);
catch ME
    pass(j) = 0;
end
j = j+1;
end
end
end






wave_types = {{'Tsunami', 9, 1, .88, [0, 0, 1], 20}};
dimensions = {{3,{'Torus',.5,1}}};
comp_type_freqs = {{'known', [], []},{'basic', 10, []},{'window one', 10, 5},{'window all', 10, 5}};
solver_types = {'basic'};

j = 1;

for i1 = 1:length(dimensions)
for i0 = 1:length(solver_types)
for i2 = 1:length(comp_type_freqs)
for i3 = 1:length(wave_types)

inobj.wave_type = wave_types{i3};
dimension = dimensions{i1}{1};
inobj.bnd_type = dimensions{i1}{2} ;
inobj.comp_type_freq = comp_type_freqs{i2};
inobj.solver_type = solver_types{i0};

s = setup_boundary(inobj.bnd_type); % may require CHEBFUN ON PATH.

inobj.s = s;
inobj.gam = s.Z;

[~, Uinc, ~, ~, ~, ~] = setup_incident(inobj.wave_type, inobj.comp_type_freq, dimension, s);

B = @(z,w) -Uinc{1}(z,w);

inobj.NystM = 75;
M = 10;


s = quadr_doubleptr_patch2(s, [M,M]);

inobj.s = s;

name{j} = [num2str(dimension),'D with ', wave_types{i3}{1} ,' wave using ', solver_types{i0} ,' solver and ', comp_type_freqs{i2}{1} ,' Fourier transform without paneling'];
pass(j) = 1;
try
    HH_density(2, B,s, inobj.solver_type);
catch ME
    pass(j) = 0;
end
j = j+1;
end
end
end
end





end
