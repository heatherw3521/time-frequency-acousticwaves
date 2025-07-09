function [pass, name] = test_plot_movie_freq_3D(obj)
% test functionality of plot_slice_freq
% this will not catch whether the plot is "wrong" in some way, just
% checks for bugs that break the code and throw an error. Heavely relies on
% plot_slive_freq_3D, less extensive testing is therefore
% done.

style = 'real';
j = 1;
pass(j) = 1;
name{j} = '3D freq movie';
file_name = 'freq_test';
try 
    plot_movie_freq_3D(obj, 20, style, file_name);
    outgiffile = fullfile(pwd, 'BroadbandHH Results/VTK/freq_test_*');
    delete(outgiffile)
catch ME
    pass(j) = 0;
end