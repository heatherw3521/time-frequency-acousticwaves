clear
%% HELP FUNCTIONS
fprintf('\n \n Help Functions.\n \n');
[p,n] = test_Clen_Curt();

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: Clen_Curt, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: Clen_Curt, %s\n', n{i1});   % Red text
    end
end


[p,n] = test_freq2time_fastsinc();

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: freq2time_fastsinc, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: freq2time_fastsinc, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_Window_ft();

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: Window_ft, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: Window_ft, %s\n', n{i1});   % Red text
    end
end

%% SETUP

fprintf('\n\n Setup Functions.\n\n');


[p,n] = test_autosample_domain();

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: autosample_domain, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: autosample_domain, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_setup_boundary();

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: setup_boundary, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: setup_boundary, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_setup_incident();

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: setup_incident, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: setup_incident, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_setup_tsunami();

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: setup_tsunami, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: setup_tsunami, %s\n', n{i1});   % Red text
    end
end

%% SOLVERS

fprintf('\n\n Solver functions:\n NOTE: some depend on eachother.\n\n');

[p,n] = test_ns_setupquad();

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: ns_setupquad, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: ns_setupquad, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_HH_density();

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: HH_density, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: HH_density, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_solve_HH_density();

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: solve_HH_density, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: solve_HH_density, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_eval_sol_freqs();

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: eval_sol_freqs, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: eval_sol_freqs, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_solve_HH_eqn();

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: solve_HH_eqn, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: solve_HH_eqn, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_eval_sol_time();

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: eval_sol_time, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: eval_sol_time, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_solve_HH_to_wave();

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: solve_HH_to_wave, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: solve_HH_to_wave, %s\n', n{i1});   % Red text
    end
end

%% 2D PLOTTING

fprintf('\n\n Plotting Functions.\n\n ');


[p,n] = test_plot_scatterer3();

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: plot_scatterer3, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: plot_scatterer3, %s\n', n{i1});   % Red text
    end
end


load("tests/functionality_tests/plotting_functions/example sol 2D.mat")

[p,n] = test_plot_slice_freq(sol);

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: plot_slice_freq, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: plot_slice_freq, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_plot_movie_freq(sol);

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: plot_movie_freq, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: plot_movie_freq, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_plot_slice_time(sol);

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: plot_slice_time, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: plot_slice_time, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_plot_movie_time(sol);

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: plot_movie_time, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: plot_movie_time, %s\n', n{i1});   % Red text
    end
end

%% 3D PLOTTING
load("tests/functionality_tests/plotting_functions/example sol 3D.mat")

[p,n] = test_plot_slice_freq_3D(sol);

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: plot_slice_freq_3D, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: plot_slice_freq_3D, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_plot_movie_freq_3D(sol);

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: plot_movie_freq_3D, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: plot_movie_freq_3D, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_plot_slice_time_3D(sol);

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: plot_slice_time_3D, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: plot_slice_time_3D, %s\n', n{i1});   % Red text
    end
end

[p,n] = test_plot_movie_time_3D(sol);

for i1 = 1:length(p)
    if p(i1)
        fprintf('Succes: plot_movie_time_3D, %s\n', n{i1}); % Green text
    else
        fprintf('Failed: plot_movie_time_3D, %s\n', n{i1});   % Red text
    end
end



