function [pass, name] = test_solve_HH_eqn()
% test functionality of solve_HH_eqn
% this will NOT catch whether the auto sampler is "wrong" in some way, just
% checks for bugs that break the code and throw an error.

field_types = {'scatter','incident', 'total'};
solver_types = {'basic'};
eval_type_freqs = {'basic', 'fmm'};
dimensions = {{2,{'Kite'}},{3,{'Torus',.5,1}}};
eval_type_times = {'freq only','fastsinc', 'GLquad', 'complexify'};


inobj.wave_type = {'Tsunami', 9, 1, .88, 1+.75*1i, 20};

inobj.delt = .1;
inobj.szIFT = 10;
inobj.NystM = 10;
inobj.bandlimit = [1,2];

inobj.comp_type_freq = {'known', [], []};


% No Paneling test
j = 1;
for i3 = 1:length(dimensions)
for i2 = 1:length(eval_type_freqs)
for i1 = 1:length(field_types)
    for i4 = 1:length(solver_types)
        for i5 = 1:length(eval_type_times)

    if i3 == 2 && i2 == 2
        continue
    end

    inobj.field_type = field_types{i1};
inobj.eval_type_freq = eval_type_freqs{i2};
inobj.solver_type = solver_types{i4};
dimension = dimensions{i3}{1};
inobj.bnd_type = dimensions{i3}{2} ;
inobj.eval_type_time = eval_type_times{i5};
inobj.eval_freqs = linspace(1,2,10);

inobj.solver_types = {'basic'};
inobj.eval_type_freqs = {'basic', 'fmm'};

if dimension == 2
    inobj.dom_range =  {[-5,5],[-5,5]};
    inobj.wave_type = {'Tsunami', 9, 1, .88, 1+.75*1i, 20};
else
    inobj.dom_range =  {[-5,5],[-5,5],[-5,5]};
    inobj.wave_type = {'Tsunami', 9, 1, .88, [0, 0, 1], 20};
end

inobj.comp_type_freq = {'known', [], []};

s = setup_boundary(inobj.bnd_type); % may require CHEBFUN ON PATH.
inobj.s = s;

[~, Uinc, ~, ~, ~, ~] = setup_incident(inobj.wave_type, inobj.comp_type_freq, dimension, s);

inobj.incF = Uinc ;

n = 10; 
inobj.autosample = {1, n}; % 1 = call autosampler, 0 = do not call
autosample = inobj.autosample{1};

if autosample
    rnge = inobj.dom_range; inobj.dom_range = rnge;
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


name{j} = [num2str(dimensions{i3}{1}),'D in ',field_types{i1},' field using ',eval_type_freqs{i2},' frequency evaluation, ',solver_types{i4},' solver and ', eval_type_times{i5}, ' inverse Fourier transform method without paneling'];
pass(j) = 1;

try
    solve_HH_eqn(inobj,zzk,0);
catch ME
    pass(j) = 0;
end
j = j+1;
end
end
end
end
end




% Paneling test

inobj.NystM = {10,2};


i3 = 1;
    for i2 = 1:length(eval_type_freqs)
    for i1 = 1:length(field_types)
        for i4 = 1:length(solver_types)
            for i5 = 1:length(eval_type_times)
    
        if i3 == 2 && i2 == 2
            continue
        end
    
        inobj.field_type = field_types{i1};
    inobj.eval_type_freq = eval_type_freqs{i2};
    inobj.solver_type = solver_types{i4};
    dimension = dimensions{i3}{1};
    inobj.bnd_type = dimensions{i3}{2} ;
    inobj.eval_type_time = eval_type_times{i5};
    inobj.eval_freqs = linspace(1,2,10);
    
    inobj.solver_types = {'basic'};
    inobj.eval_type_freqs = {'basic', 'fmm'};
    
    if dimension == 2
        inobj.dom_range =  {[-5,5],[-5,5]};
        inobj.wave_type = {'Tsunami', 9, 1, .88, 1+.75*1i, 20};
    else
        inobj.dom_range =  {[-5,5],[-5,5],[-5,5]};
        inobj.wave_type = {'Tsunami', 9, 1, .88, [0, 0, 1], 20};
    end
    
    inobj.comp_type_freq = {'known', [], []};
    
    s = setup_boundary(inobj.bnd_type); % may require CHEBFUN ON PATH.
    inobj.s = s;

    [~, Uinc, ~, ~, ~, ~] = setup_incident(inobj.wave_type, inobj.comp_type_freq, dimension,s);
    
    inobj.incF = Uinc ;
    
    n = 10; 
    inobj.autosample = {1, n}; % 1 = call autosampler, 0 = do not call
    autosample = inobj.autosample{1};
    
    if autosample
        rnge = inobj.dom_range; inobj.dom_range = rnge;
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
    
    
    name{j} = [num2str(dimensions{i3}{1}),'D in ',field_types{i1},' field using ',eval_type_freqs{i2},' frequency evaluation, ',solver_types{i4},' solver and ', eval_type_times{i5}, ' inverse Fourier transform method with paneling'];
    pass(j) = 1;
    
    try
        solve_HH_eqn(inobj,zzk,0);
    catch ME
        pass(j) = 0;
    end
    j = j+1;
    end
    end
    end
end

% PLotting test

inobj.field_type = field_types{i1};
inobj.eval_type_freq = eval_type_freqs{i2};
inobj.solver_type = solver_types{i4};
dimension = dimensions{i3}{1};
inobj.bnd_type = dimensions{i3}{2} ;
inobj.eval_type_time = eval_type_times{i5};
inobj.eval_freqs = linspace(1,2,10);

inobj.solver_types = {'basic'};
inobj.eval_type_freqs = {'basic', 'fmm'};

if dimension == 2
    inobj.dom_range =  {[-5,5],[-5,5]};
    inobj.wave_type = {'Tsunami', 9, 1, .88, 1+.75*1i, 20};
else
    inobj.dom_range =  {[-5,5],[-5,5],[-5,5]};
    inobj.wave_type = {'Tsunami', 9, 1, .88, [0, 0, 1], 20};
end

inobj.comp_type_freq = {'known', [], []};

s = setup_boundary(inobj.bnd_type); % may require CHEBFUN ON PATH.
inobj.s = s;

[~, Uinc, ~, ~, ~, ~] = setup_incident(inobj.wave_type, inobj.comp_type_freq, dimension, s);

inobj.incF = Uinc ;

n = 10; 
inobj.autosample = {1, n}; % 1 = call autosampler, 0 = do not call
autosample = inobj.autosample{1};

if autosample
    rnge = inobj.dom_range; inobj.dom_range = rnge;
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


name{j} = ['panel plotting'];
pass(j) = 1;

try
    solve_HH_eqn(inobj,zzk,1);
catch ME
    pass(j) = 0;
end
j = j+1;
clf; close()
