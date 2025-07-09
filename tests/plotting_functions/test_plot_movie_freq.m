function [pass, name] = test_plot_movie_freq(obj)
% test functionality of plot_movie_freq
% this will not catch whether the plot is "wrong" in some way, just
% checks for bugs that break the code and throw an error. Heavely relies on
% plot_slive_freq & plot_scatter3, less extensive testing is therefore
% done.

j = 1;
pass(j) = 1;
name{j} = 'freq movie';
plottype = {'mesh', 2}; style = 'real'; rng = []; clr = []; 
try
    plot_movie_freq(obj, 5, style, plottype, rng, clr, 'test_movie_freq',1)
    outgiffile = fullfile(pwd, 'BroadbandHH Results/Gifs/test_movie_freq.avi');
    delete(outgiffile)
    clf
    close all
catch ME
    pass(j) = 0;
end
