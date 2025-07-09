function [pass, name] = test_plot_slice_freq_3D(obj)
% test functionality of plot_slice_freq_3D
% this will not catch whether the plot is "wrong" in some way, just
% checks for bugs that break the code and throw an error

k = 8; style = 'real';
j = 1;
pass(j) = 1;
name{j} = '3D freq slice';
file_name = 'freq_test';
try 
    plot_slice_freq_3D(obj, k, style, file_name);
    outgiffile = fullfile(pwd, 'BroadbandHH Results/VTK/freq_test_8.vtk');
    delete(outgiffile)
catch ME
    pass(j) = 0;
end