function [pass, name] = test_plot_movie_time(obj)
% test functionality of plot_slice_freq
% this will not catch whether the plot is "wrong" in some way, just
% checks for bugs that break the code and throw an error. Heavely relies on
% plot_slive_time & plot_scatter3, less extensive testing is therefore
% done.

j = 1;
pass(j) = 1;
name{j} = 'time movie';
plottype = {'mesh', 2}; style = 'real'; rng = []; clr = []; 
try
    plot_movie_time(obj, 5, style, plottype, rng, clr, 'test_movie_time',1)
    outgiffile = fullfile(pwd, 'BroadbandHH Results/Gifs/test_movie_time.mp4');
    delete(outgiffile)
    clf
    close all
catch ME
    pass(j) = 0;
end