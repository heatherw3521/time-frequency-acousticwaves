function  plot_slice_time(obj, k,style,plottype, rng, clr, im_name, save_im)
 %
 % plot the solution slice at time indexed by k
 % over all points in stored in obj as its domain. 
 %
 % style = 'abs', 'real', 'imag'
 %
 % rng = [xlim, ylim]; 
 % clr = colormap. 
 % plottype = {mesh or surf, 2 = view(2), 3 = view(3)}; 
 %%
 

vals = obj.timevals; vals = vals{1}; % pull out of cell
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
    rngz = []; 
end
if length(rng)==1
    rngz = rng{:};
    rng = obj.dom_range; 
    rngx = rng{1};
    rngy = rng{2};
elseif length(rng)==2
    rngx = rng{1}; 
    rngy = rng{2}; 
    rngz = []; 
else
    rngx = rng{1}; 
    rngy = rng{2}; 
    rngz = rng{3};
end

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
if ~isempty(rngz)
    zlim(rngz)
end
set(gcf,'color','w')
t = obj.times;
title(['Solution at time $t =~$',num2str(t(k),'%.3f')],'Interpreter','latex','FontSize',20)
set(gca, 'fontsize', 18)
axis equal
axis off
axis tight

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
