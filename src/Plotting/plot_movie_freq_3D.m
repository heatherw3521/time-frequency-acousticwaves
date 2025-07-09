function plot_movie_freq_3D(obj, sp, style, file_name)
vals = obj.freqvals; vals = vals{1}; % pull out of cell
n = size(vals,2);
for k = sp:n
    plot_slice_freq_3D(obj, k, style, file_name)
end
disp(['VTK files saved at ',fullfile(pwd, ['BroadbandHH Results/VTK/',file_name]),'.'])
end