function [pass, name] = test_plot_slice_freq(obj)
% test functionality of plot_slice_freq
% this will not catch whether the plot is "wrong" in some way, just
% checks for bugs that break the code and throw an error

%%
% mesh plots  #1
k = 5; 
j = 1; 
pass(j) = 1; 
name{j} = 'mesh view 2';
plottype = {'mesh', 2}; style = 'real'; rng = []; clr = []; 
try
plot_slice_freq(obj, k,style,plottype, rng, clr,"a",0)
catch ME
    pass(j) = 0;
end
j = j+1;
%          #2
pass(j) = 1;
name{j} = 'mess view 3';
plottype = {'mesh', 3};
try
plot_slice_freq(obj, k,style,plottype, rng, clr,"a",0)
catch ME
    pass(j) = 0;
end
j = j+1; 
%%
% surf plots:  #3 
pass(j) = 1; 
name{j} = "surf view 2";
plottype = {'mesh', 2}; style = 'real'; rng = []; clr = []; 
try
plot_slice_freq(obj, k,style,plottype, rng, clr,"a",0)
catch ME
    pass(j) = 0;
end
j = j+1;
%            #4
pass(j) = 1;
name{j} = 'surf view 3';
plottype = {'mesh', 3};
try
plot_slice_freq(obj, k,style,plottype, rng, clr,"a",0)
catch ME
    pass(j) = 0;
end
j = j+1; 
%%
% pass in rng:  #5
pass(j) = 1; 
name{j} = 'new range'; 
rng = {[-2, 2],[ -3, 3]}; 
try
plot_slice_freq(obj, k,style,plottype, rng, clr,"a",0)
catch ME
    pass(j) = 0;
end
j = j+1; 
%%
% pass in colormap (requires cmocean on path)   #6
pass(j) = 1; 
name{j} = 'colormap parula string';
rng = []; 
clr = "parula"; 
try
plot_slice_freq(obj, k,style,plottype, rng, clr,"a",0)
catch ME
    pass(j) = 0;
end
j = j+1; 
%                         #7
pass(j) = 1; 
name{j} = 'colormap matrix'; 
rng = []; 
clr = cmocean('curl'); 
try
plot_slice_freq(obj, k,style,plottype, rng, clr,"a",0)
catch ME
    pass(j) = 0;
end
clf; close();
end



