function [pass, name] = test_setup_boundary()
% test functionality of plot_setup_incident
% this will NOT catch whether the boundary setup is "wrong" in some way, just
% checks for bugs that break the code and throw an err

Boundary_names_2D = {{'Kite'},{'Star',.3,7},{'Circles',[0 + 0i, 4*(-0.5+1*1i), 5*(1+.5*1i)^(-1)],[2,1,1]},{'C Curve',3,2.5,1,exp(-1i*3*pi/4)},{'Whisper Gallery',3,6/5*3,1.1},{'Crescent Gallery',.1,.9,4},{'Keyhole',1,3,.5,0}, {'HouseHolder'}, {'HouseHolder_sky'}, {'Teardrop', 1},{'Radiator',1.5,5,.5,.8,.2,3,pi/2},{'Spikes',3},{'Von Koch Snowflake', 3},{'Perforated Wall',5,.1,pi/4,6}};
j = 1;
for i1 = 1:length(Boundary_names_2D)
    name{i1} = Boundary_names_2D{i1}{1};
    pass(j) = 1;
    try
        s = setup_boundary(Boundary_names_2D{i1});
        s.Z(1);
    catch ME
        pass(j) = 0;
    end
    j = j+1;
end


Boundary_names_3D = {{'Torus',.5,1},{'Wobbly torus',3,2,.2},{'C Torus',5,2.5,.5,exp(1i*pi/4),1}};

n = 10;
for i1 = 1:length(Boundary_names_3D)
    pass(j) = 1;    
    name{j} = Boundary_names_3D{i1}{1};
    try
        s = setup_boundary(Boundary_names_3D{i1});
        s.Z(1,2);
    catch ME
        pass(j) = 0;
    end
    j = j+1;
end
