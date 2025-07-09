function [pass, name] = test_setup_incident()
% test functionality of plot_setup_incident
% this will NOT catch whether the auto sampler is "wrong" in some way, just
% checks for bugs that break the code and throw an error

j = 1;

%% Tsunami
wave_type = {'Tsunami', 9, 1, .88, 1+.75*1i, 20};
comp_type_freq = {'known'};
dimension = 2;

[pass(j), name{j}] = try_catch_gen(wave_type, comp_type_freq, dimension);
j = j + 1;


wave_type = {'Tsunami', 9, 1, .88, 1+.75*1i, 20};
comp_type_freq = {'basic',50};
dimension = 2;

[pass(j), name{j}] = try_catch_gen(wave_type, comp_type_freq, dimension);
j = j + 1;


wave_type = {'Tsunami', 9, 1, .88, 1+.75*1i, 20};
comp_type_freq = {'window all', 50, 10};
dimension = 2;

[pass(j), name{j}] = try_catch_gen(wave_type, comp_type_freq, dimension);
j = j + 1;


wave_type = {'Tsunami', 9, 1, .88, 1+.75*1i, 20};
comp_type_freq = {'window one', 50, 10};
dimension = 2;

[pass(j), name{j}] = try_catch_gen(wave_type, comp_type_freq, dimension);
j = j + 1;


wave_type = {'Tsunami', 9, 1, .88, [0, 0, 1], 20};
comp_type_freq = {'known'};
dimension = 3;

[pass(j), name{j}] = try_catch_gen(wave_type, comp_type_freq, dimension);
j = j + 1;


wave_type = {'Tsunami', 9, 1, .88, [0, 0, 1], 20};
comp_type_freq = {'basic',50};
dimension = 3;

[pass(j), name{j}] = try_catch_gen(wave_type, comp_type_freq, dimension);
j = j + 1;


wave_type = {'Tsunami', 9, 1, .88, [0, 0, 1], 20};
comp_type_freq = {'window all', 50, 10};
dimension = 3;

[pass(j), name{j}] = try_catch_gen(wave_type, comp_type_freq, dimension);
j = j + 1;


wave_type = {'Tsunami', 9, 1, .88, [0, 0, 1], 20};
comp_type_freq = {'window one', 50, 10};
dimension = 3;

[pass(j), name{j}] = try_catch_gen(wave_type, comp_type_freq, dimension);
j = j + 1;

end

function [pass, name] = try_catch_gen(wave_type, comp_type_freq, dimension)
pass = 1;
name = [num2str(dimension),'D ',wave_type{1}, ' (',comp_type_freq{1},')'];

if dimension == 2
    s = setup_boundary({'Kite'});
else
    s = setup_boundary({'Torus',.5,1});
end

try
    [uinc, Uinc, Uinc_slow, ~, ~,~] = setup_incident(wave_type, comp_type_freq, dimension, s);

    if dimension == 2
        uinc{1}( [.21+ .3i; .23+ .5i],2);
        Uinc{1}( [.21+ .3i; .23+ .5i],.2i);
        if ~isempty(Uinc_slow)
            Uinc_slow{3}( [.21+ .3i; .23+ .5i],.2i);
        end
    else
        uinc{1}( [.21+ .3i,1,1; .23+ .5i,1,1;.23 - .5i,1,1],2);
        Uinc{1}( [.21+ .3i,1,1; .23+ .5i,1,1;.23 - .5i,1,1],.2i);
    end
catch ME
    pass = 0;
end
end