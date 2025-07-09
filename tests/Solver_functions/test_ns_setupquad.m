function [pass, name] = test_ns_setupquad()

    Boundary_names = {{'Kite'},{'Star',.3,7},{'Circles',[0 + 0i, 4*(-0.5+1*1i), 5*(1+.5*1i)^(-1)],[2,1,1]},{'C Curve',3,2.5,1,exp(-1i*3*pi/4)},{'Whisper Gallery',3,6/5*3,1.1},{'Crescent Gallery',.1,.9,4},{'Keyhole',1,3,.5,0}, {'HouseHolder'}, {'HouseHolder_sky'}, {'Teardrop', 1},{'Radiator',1.5,5,.5,.8,.2,3,pi/2},{'Spikes',3},{'Von Koch Snowflake', 2},{'Perforated Wall',5,.1,pi/4,6}};
    j = 1;



    M = {100,0};

    plotme = 0;
    
    for i1 = 1:length(Boundary_names)
        s = setup_boundary(Boundary_names{i1});

        pass(j) = 1;

        name{j} = ['No subdivisions with no plotting over scatterer ',Boundary_names{i1}{1}];

        try
            ns_setupquad(s, M, plotme);
        catch ME
            pass(j) = 0;
        end

        j = j+1;
    end
    
    M = {300,2};
    
    plotme = 0;

    for i1 = 1:length(Boundary_names)
        s = setup_boundary(Boundary_names{i1});

        pass(j) = 1;

        name{j} = ['Subdivisions with no plotting over scatterer ',Boundary_names{i1}{1}];

        try
            ns_setupquad(s, M, plotme);
        catch ME
            pass(j) = 0;
        end

        j = j+1;
    end
    

    M = {100,0};

plotme = 1;

for i1 = 1:length(Boundary_names)
    s = setup_boundary(Boundary_names{i1});

    pass(j) = 1;

    name{j} = ['No subdivisions with plotting over scatterer ',Boundary_names{i1}{1}];

    try
        ns_setupquad(s, M, plotme);
        clf; close();
    catch ME
        pass(j) = 0;
    end

    j = j+1;
end

M = {300,2};

plotme = 1;

for i1 = 1:length(Boundary_names)
    s = setup_boundary(Boundary_names{i1});

    pass(j) = 1;

    name{j} = ['Subdivisions with plotting over scatterer ',Boundary_names{i1}{1}];

    try
        ns_setupquad(s, M, plotme);
        clf; close();
    catch ME
        pass(j) = 0;
    end

    j = j+1;
end


