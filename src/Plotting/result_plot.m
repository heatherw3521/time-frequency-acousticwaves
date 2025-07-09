function result_plot(inobj,vidobj)

%ADD TEMPLATES

dimension = 2;
if strcmp(inobj.bnd_type{1},'Torus') || strcmp(inobj.bnd_type{1},'Wobbly torus') || strcmp(inobj.bnd_type{1},'Pot') || strcmp(inobj.bnd_type{1},'C Torus')
    dimension = 3;
end

if dimension == 2

        if strcmp(vidobj.soltype,'freq')
            if strcmp(inobj.eval_type_time,'complexify')
                inobj.freqvals = {inobj.freqvals{1}{2}};
                inobj.freq = {inobj.freq{2}};
            end
            vidobj.endpoint = min(vidobj.endpoint,length(inobj.freqs));

            if vidobj.frame == 0
                plot_movie_freq(inobj, vidobj.startpoint, vidobj.valuetype, vidobj.plottype,vidobj.rng, vidobj.clr, vidobj.im_name,vidobj.saveim)
            else
                plot_slice_freq(inobj, vidobj.frame, vidobj.valuetype,vidobj.plottype, vidobj.rng, vidobj.clr, vidobj.im_name,vidobj.saveim)
            end

        else
                        vidobj.endpoint = min(vidobj.endpoint,length(inobj.times));

            if vidobj.frame == 0
                plot_movie_time(inobj, vidobj.startpoint, vidobj.endpoint, vidobj.valuetype,vidobj.plottype, vidobj.rng, vidobj.clr, vidobj.im_name,vidobj.saveim)
            else
                plot_slice_time(inobj, vidobj.frame, vidobj.valuetype,vidobj.plottype, vidobj.rng, vidobj.clr, vidobj.im_name,vidobj.saveim)
            end
            
        end

else

    if vidobj.saveim
        if strcmp(vidobj.soltype,'freq')
            if strcmp(inobj.eval_type_time,'complexify')
                inobj.freqvals = {inobj.freqvals{1}{2}};
                inobj.freq = {inobj.freq{2}};
            end

            if vidobj.frame == 0
                plot_movie_freq_3D(inobj, vidobj.startpoint, vidobj.valuetype, vidobj.im_name)
            else
                plot_slice_freq_3D(inobj, vidobj.frame, vidobj.valuetype, vidobj.im_name)
            end

        else

            if vidobj.frame == 0
                plot_movie_time_3D(inobj, vidobj.startpoint, vidobj.valuetype, vidobj.im_name)
            else
                plot_slice_time_3D(inobj, vidobj.frame, vidobj.valuetype, vidobj.im_name)
            end

        end
    end

end

end