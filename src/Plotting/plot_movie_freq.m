function plot_movie_freq(obj, sp, style, plottype, rng, clr, Gif_name, save_vid)
%
% Create mp4 video with name Gif_name of
% the solution from index sp onward
% over all points in stored in obj as its domain.
% NOTE: must be in folder "Faster Clean".
%
% style = 'abs', 'real', 'imag'
%
% rng = [xlim, ylim];
% clr = colormap.
% plottype = {mesh or surf, 2 = view(2), 3 = view(3)};

vals = obj.freqvals; vals = vals{1}; % pull out of cell

if length(rng) == 2 || length(rng) == 0
    if strcmpi(style, 'abs')
        maxvals = max(max(abs(vals)));
        minvals = min(min(abs(vals)));
    elseif strcmpi(style, 'real')
        maxvals = max(max(real(vals)));
        minvals = min(min(real(vals)));
    elseif strcmpi(style, 'imag')
        maxvals = max(max(imag(vals)));
        minvals = min(min(imag(vals)));
    else
        error('real, imag, or abs plot?')
    end
    zrng = [minvals, maxvals];
elseif length(rng) == 1
    zrng = rng{1};
    rng = [];
else
    zrng = rng{3};
end

t = obj.freq;
clf
close all


figure(365); hold on
for k = sp:length(t)
    plot_slice_freq(obj, k, style, plottype, rng, clr, Gif_name, 0)

    zlim(zrng);
    caxis(zrng);

    plot_scatterer3(obj.gam, .5,1) % This can be commented out if you don't want the scatterer
    title(['Solution at frequency $k =~$',num2str(t(k),'%.3f')],'Interpreter','latex','FontSize',20)

    if save_vid
        Image = getframe(gcf);

        outputDir = fullfile(pwd, 'BroadbandHH Results','Gifs');
        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        end

        outputFile = fullfile(outputDir, ['im_', num2str(k - sp), '.png']);

        imwrite(Image.cdata, outputFile,'png');
    end
    clf;
end
if save_vid
    Gif_writing(Gif_name);
    disp(['Movie saved at ',outputDir,'.'])
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%subfunctions%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GIF maker

function [] = Gif_writing(gifFile)
outputDir = fullfile(pwd, 'BroadbandHH Results/Gifs');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
outgifFile = fullfile(pwd, ['BroadbandHH Results/Gifs/',gifFile]);
i = 1;
while isfile(fullfile(outputDir, ['im_', num2str(i - 1), '.png']))
    images{i} = imread(fullfile(outputDir, ['im_', num2str(i - 1), '.png']));
    delete(fullfile(outputDir, ['im_', num2str(i - 1), '.png']))
    i = i + 1;
end

n = i-1;
writerObj = VideoWriter(outgifFile,'Motion JPEG AVI');
%writerObj = VideoWriter(gifFile,'Uncompressed AVI');
%writerObj.FrameRate = 30;
writerObj.FrameRate = 50;
% set the seconds per image
secsPerImage = ones(n) ;
% open the video writer
open(writerObj);
% write the frames to the video
for u=1:n
    % convert the image to a frame
    frame = im2frame(images{u});
    for v=1:secsPerImage(u)
        writeVideo(writerObj, frame);
    end
end
% close the writer object
close(writerObj);
end