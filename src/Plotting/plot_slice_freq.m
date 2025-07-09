function  plot_slice_freq(obj, k,style,plottype, rng, clr, im_name, save_im)
 %
 % plot the solution slice at freq indexed by k
 % over all points in stored in obj as its domain. 
 %
 % style = 'abs', 'real', 'imag'
 %
 % rng = [xlim, ylim]; 
 % clr = colormap. 
 % plottype = {mesh or surf, 2 = view(2), 3 = view(3)}; 
 %%


vals = obj.freqvals; vals = vals{1}; % pull out of cell
t = obj.freq;
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
% check for color
if isempty(clr)
    clr = "cool"; %set to default MATLAB color
end

if isempty(rng)
    rng = obj.dom_range; 
end
rngx = rng{1}; 
rngy = rng{2}; 

if isempty(plottype)
    issurf = 1; 
    viewtype = 2; 
else
    layer = plottype{1}; 
    viewtype = plottype{2}; 
    if strcmpi(layer, 'mesh')
        issurf = 0; 
    else 
        issurf = 1; 
    end
end

% unpack domain: 
domx = obj.domaingrid{1};
domy = obj.domaingrid{2};  
keepfilter = obj.domaingrid{4};
sz = size(domx);

pltvals = zeros(sz); 
pltvals(keepfilter) = vals; pltvals(~keepfilter) = nan;
pltvals = reshape(pltvals, sz); 

if issurf
    surf(domx, domy, pltvals)
else
    mesh(domx, domy, pltvals)
end
view(viewtype)
shading flat
hold on
colormap(clr)
colorbar
xlim(rngx); ylim(rngy);
set(gcf,'color','w')
title(['Solution at frequency $k =~$',num2str(t(k),'%.3f')],'Interpreter','latex','FontSize',20)
set(gca, 'fontsize', 18)

if save_im
    Image = getframe(gcf);
    outputDir = fullfile(pwd, 'BroadbandHH Results','Images');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    outputFile = fullfile(outputDir, [im_name,'_', num2str(k), '.png']);
    imwrite(Image.cdata, outputFile,'png');
end

hold off

end
