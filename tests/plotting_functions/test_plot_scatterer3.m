function [pass, name] = test_plot_scatterer3()
% test functionality of plot_plot_scatterer3
% this will not catch whether the plot is "wrong" in some way, just
% checks for bugs that break the code and throw an error

%TODO: 3D plots

Boundary_names_2D = {{'Kite'},{'Star',.3,7},{'Circles',[0 + 0i, 4*(-0.5+1*1i), 5*(1+.5*1i)^(-1)],[2,1,1]},{'C Curve',3,2.5,1,exp(-1i*3*pi/4)},{'Whisper Gallery',3,6/5*3,1.1},{'Crescent Gallery',.1,.9,4},{'Keyhole',1,3,.5,pi}, {'Teardrop', 1},{'Radiator',1.5,5,.5,.8,.2,3,pi/2},{'Spikes',3},{'Von Koch Snowflake', 3},{'Perforated Wall',5,.1,pi/4,6}};
j = 1;
for i1 = 1:length(Boundary_names_2D)
    name{j} = Boundary_names_2D{i1}{1};
    pass(j) = 1;
    s = setup_boundary(Boundary_names_2D{i1});
    try
        plot_scatterer3(s.Z, 0, 0); clf; close all
    catch ME
        pass(j) = 0;
    end
    j = j+1;
end

for i1 = 1:length(Boundary_names_2D)
    name{j} = ['filled ', Boundary_names_2D{i1}{1}];
    pass(j) = 1;
    s = setup_boundary(Boundary_names_2D{i1});
    try
        plot_scatterer3(s.Z, 0, 1); clf; close all
    catch ME
        pass(j) = 0;
    end
    j = j+1;
end