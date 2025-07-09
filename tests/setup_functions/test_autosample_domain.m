function [pass, name] = test_autosample_domain()
% test functionality of plot_autosample_domain
% this will NOT catch whether the auto sampler is "wrong" in some way, just
% checks for bugs that break the code and throw an error

j = 1;
Boundary_names_2D = {{'Kite'},{'Star',.3,7},{'Circles',[0 + 0i, 4*(-0.5+1*1i), 5*(1+.5*1i)^(-1)],[2,1,1]},{'C Curve',3,2.5,1,exp(-1i*3*pi/4)},{'Whisper Gallery',3,6/5*3,1.1},{'Crescent Gallery',.1,.9,4},{'Keyhole',1,3,.5,0}, {'Teardrop', 1},{'Radiator',1.5,5,.5,.8,.2,3,pi},{'Spikes',3},{'Von Koch Snowflake', 3},{'Perforated Wall',5,.1,pi/4,6}};
rnge = {[-5,6],[-6,4.5]};
n = 10;
for i1 = 1:length(Boundary_names_2D)
    pass(j) = 1;    
    name{j} = Boundary_names_2D{i1}{1};
    try
        [~,~,~,~] = autosample_domain(n, rnge, Boundary_names_2D{i1}, 0);
    catch ME
        pass(j) = 0;
    end
    j = j+1;
end

Boundary_names_3D = {{'Torus',.5,1},{'Wobbly torus',3,2,.2},{'C Torus',5,2.5,.5,exp(1i*pi/4),1}};
rnge = {[-5,6],[-6,4.5],[-5,5]};
n = 10;
for i1 = 1:length(Boundary_names_3D)
    pass(j) = 1;    
    name{j} = Boundary_names_3D{i1}{1};
    try
        [~,~,~,~] = autosample_domain(n, rnge, Boundary_names_3D{i1}, 0);
    catch ME
        pass(j) = 0;
    end
    j = j+1;
end

