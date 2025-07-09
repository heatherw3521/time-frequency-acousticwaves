function plot_slice_freq_3D(obj, k, style, file_name)

vals = obj.freqvals; vals = vals{1}; % pull out of cell
vals = vals(:,k);
if strcmpi(style, 'abs')
    vals = abs(vals); 
elseif strcmpi(style, 'real')
    vals = real(vals);
elseif strcmpi(style, 'imag')
    vals = imag(vals);
else
    error('real, imag, or abs plot?')
end

rng = obj.dom_range;
rngx = rng{1};
rngy = rng{2};
rngz = rng{3};

domaingrid = obj.domaingrid;
xx = domaingrid{1};
keep = domaingrid{4};
n = round(length(xx)^(1/3));

[xx,yy,zz] = meshgrid(linspace(rngx(1),rngx(2),n),linspace(rngy(1),rngy(2),n),linspace(rngz(1),rngz(2),n));

pltvals = zeros(n,1);
pltvals(keep) = vals;
pltvals = reshape(pltvals,n,n,n);

outputDir = fullfile(pwd, 'BroadbandHH Results','VTK');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

vtk_name = fullfile(outputDir, [file_name,'_',num2str(k),'.vtk']);

vtkwrite(vtk_name,'structured_grid',xx,yy,zz,'scalars','Amplitude',pltvals);