function [pass, name] = test_eval_sol_freqs()
% test functionality of eval_sol_freqs
% this will NOT catch whether the auto sampler is "wrong" in some way, just
% checks for bugs that break the code and throw an error

field_types = {'scatter','incident', 'total'};
eval_type_freqs = {'basic', 'fmm'};
dimensions = {{2,{'Kite'}},{3,{'Torus',.5,1}}};

j = 1;
for i3 = 1:length(dimensions)
for i2 = 1:length(eval_type_freqs)
for i1 = 1:length(field_types)
    clearvars -except field_types eval_type_freqs dimensions i1 i2 i3 j pass name

    if i3 == 2 && i2 == 2
        continue
    end

    field_type = field_types{i1};
eval_type_freq = eval_type_freqs{i2};
dimension = dimensions{i3}{1};
inobj.bnd_type = dimensions{i3}{2} ;

if dimension == 2
    inobj.dom_range =  {[-5,5],[-5,5]};
    inobj.wave_type = {'Tsunami', 9, 1, .88, 1+.75*1i, 20};
else
    inobj.dom_range =  {[-5,5],[-5,5],[-5,5]};
    inobj.wave_type = {'Tsunami', 9, 1, .88, [0, 0, 1], 20};
end

% Incident wave
inobj.comp_type_freq = {'known', [], []};

% obtaining s
s = setup_boundary(inobj.bnd_type); % may require CHEBFUN ON PATH.

[~, Uinc, ~, ~, ~, ~] = setup_incident(inobj.wave_type, inobj.comp_type_freq, dimension, s);

M = 15;

% Solver type
inobj.solver_type = 'basic';
ww = linspace(1,2,10);
sol_HH_dens = solve_HH_density(Uinc{1}, s,  [1,2], ww, M, inobj.solver_type);

n = 10; 
inobj.autosample = {1, n}; % 1 = call autosampler, 0 = do not call
autosample = inobj.autosample{1};

if autosample
    rnge = inobj.dom_range; sol.dom_range = rnge;
    [zz, xx, yy, keep] = autosample_domain(n, rnge, inobj.bnd_type, 0);
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


name{j} = [num2str(dimensions{i3}{1}),'D in ',field_types{i1},' field using ',eval_type_freqs{i2},' frequency evaluation without paneling'];
pass(j) = 1;
try
    eval_sol_freqs(sol_HH_dens, zzk, field_type, eval_type_freq, ww);
catch
    pass(j) = 0;
end
j = j+1;
end
end
end

% With panel test
i3 = 1;
for i2 = 1:length(eval_type_freqs)
for i1 = 1:length(field_types)
    clearvars -except field_types eval_type_freqs dimensions i1 i2 i3 j pass name

    if i3 == 2 && i2 == 2
        continue
    end

    field_type = field_types{i1};
eval_type_freq = eval_type_freqs{i2};
dimension = dimensions{i3}{1};
inobj.bnd_type = dimensions{i3}{2} ;

if dimension == 2
    inobj.dom_range =  {[-5,5],[-5,5]};
    inobj.wave_type = {'Tsunami', 9, 1, .88, 1+.75*1i, 20};
else
    inobj.dom_range =  {[-5,5],[-5,5],[-5,5]};
    inobj.wave_type = {'Tsunami', 9, 1, .88, [0, 0, 1], 20};
end

% Incident wave

inobj.wave_type = {'Tsunami', 9, 1, .88, 1+.75*1i, 20};
inobj.comp_type_freq = {'known', [], []};

% obtaining s
s = setup_boundary(inobj.bnd_type); % may require CHEBFUN ON PATH.

[~, Uinc, ~, ~, ~, ~] = setup_incident(inobj.wave_type, inobj.comp_type_freq, dimension, s);

M = {15,2};

% Solver type
inobj.solver_type = 'basic';
ww = linspace(1,2,10);
sol_HH_dens = solve_HH_density(Uinc{1}, s,  [1,2], ww, M, inobj.solver_type,0);

n = 10; 
inobj.autosample = {1, n}; % 1 = call autosampler, 0 = do not call
autosample = inobj.autosample{1};

if autosample
    rnge = inobj.dom_range; sol.dom_range = rnge;
    [zz, xx, yy, keep] = autosample_domain(n, rnge, inobj.bnd_type, 0);
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


name{j} = [num2str(dimensions{i3}{1}),'D in ',field_types{i1},' field using ',eval_type_freqs{i2},' frequency evaluation with paneling'];
pass(j) = 1;
try
    eval_sol_freqs(sol_HH_dens, zzk, field_type, eval_type_freq, ww);
catch
    pass(j) = 0;
end
j = j+1;
end
end