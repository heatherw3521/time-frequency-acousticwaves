function s = setup_boundary(bnd_type)

if isa(bnd_type{1}, 'chebfun')
    gam = bnd_type{1};
    dgam = diff(gam);
    dgam2 = diff(dgam);
    s.Z = @(t) gam(t);
    s.Zp = @(t) dgam(t);
    s.Zpp = @(t) dgam2(t);
    return
end
switch bnd_type{1}
    case 'Kite'

        a = 2;
        gam = chebfun(@(t) cos(t) + .65*cos(a*t)-.65 + 1i*(1.5*sin(t)), [0, 2*pi]);

    case 'Star'

        a = bnd_type{2}; w = bnd_type{3};
        gam = chebfun(@(t) (2 + 2*a*cos(w*t)).*exp(1i*t), [0, 2*pi]);

    case 'Circles'

        c = bnd_type{2}; % centers
        R = bnd_type{3}; %radii
        NN = length(R);
        a = 2*pi / NN;
        t = chebfun('t',[0, a],'splitting','on');
        gam = c(1) + R(1) * exp(1i*t*NN);
        for i1 = 2:length(c)
            gam = join(gam,c(i1) + R(i1) * exp(1i*t*NN));
        end

    case 'Whisper Gallery'

        % Not parabolas, need to switch it!
        a = bnd_type{2}; %
        b = bnd_type{3};
        fat = bnd_type{4};

            c = a*fat;
            d = b/fat;
            
        %Crossing points
        x = (1/a^2+(1/c^2-1/a^2)/(1-b^2/d^2))^(-1/2);
        y = ((1/c^2-1/a^2)/(1/b^2-1/d^2))^(1/2)*abs(x);
        
        %Parametrizations
        gam1 = @(t) a * cos(t) + 1i * b *sin(t);
        gam2 = @(t) c * cos(t) + 1i * d *sin(t);
        
        %Crossing tvalues chronological
        z = [ x/a - y/b*1i, x/a + y/b*1i, - x/a + y/b*1i, -x/a - y/b*1i];
        ts1 = angle(z); ts1(ts1<-pi/2) = ts1(ts1<-pi/2) + 2*pi;
        z = [x/c - y/d*1i,x/c + y/d*1i, - x/c + y/d*1i, -x/c - y/d*1i];
        ts2 = angle(z); ts2(ts2<-pi/2) = ts2(ts2<-pi/2) + 2*pi;
        
        
        aa = 2*pi/4;
        t = chebfun('t',[0, aa],'splitting','on');

        gam = a * cos(t/aa*ts1(2) + (1-t/aa)*ts1(1)) + 1i * b *sin(t/aa*ts1(2) + (1-t/aa)*ts1(1));

        gam = join(gam, c * cos(t/aa*ts2(1) + (1-t/aa)*ts2(2)) + 1i * d *sin(t/aa*ts2(1) + (1-t/aa)*ts2(2)));
        gam = join(gam, a * cos(t/aa*ts1(4) + (1-t/aa)*ts1(3)) + 1i * b *sin(t/aa*ts1(4) + (1-t/aa)*ts1(3)));
        gam = join(gam, c * cos(t/aa*ts2(3) + (1-t/aa)*ts2(4)) + 1i * d *sin(t/aa*ts2(3) + (1-t/aa)*ts2(4)));
            
    case 'C Curve'

        a = bnd_type{2}; %controls how "sharp" the corners are, bigger = sharper
        b = bnd_type{3}; %2.5 (smaller = less extreme cavity/less depth into the C)
        c = bnd_type{4}; % .1 (larger = "fatter all around")
        d = bnd_type{5}; % cavity opens to the left; choose as exp(1i*(0--2*pi)) to rotate
        t = chebfun('t', [0, 2*pi]);
        r = 3+c*tanh(a*cos(t));
        th = b*sin(t);
        gam = d*exp(1i*th).*r;

    case 'Crescent Gallery'
        r = 5; 
        th = exp(1i*pi/2);
        p = 1i; 
        a1 = bnd_type{2}; % ?
        a2 = bnd_type{3}; % ?
        d = bnd_type{4}; % Distance between crescents
    
        t = chebfun('t', [0, pi],'splitting','on');

        gam = conj(exp(2i*t)- a1/(exp(2i*t)+a2))+d/2;
        gam = join(gam+p, -conj(exp(2i*(t+pi))- a1/(exp(2i*(t+pi))+a2))-d/2);
        gam = r*th*gam;

    case 'HouseHolder'
        Houzes = [-4 + 1.54*1i,-2.8+1.54*1i,-2.8+2.03*1i,-2.35+2.03*1i,-2.35+.22*1i,-2.26+.22*1i,-2.26+2.05*1i,-1.46+2.05*1i,-1.46+.15*1i,-1.19+.15*1i,-1.19+.53*1i,-.82+.53*1i,-.49+.89*1i,-.49+.97*1i,-.31+.97*1i,-.31+.89*1i, -.01+.53*1i,1.26+.53*1i, 1.31 + 2.32*1i,1.58+2.95*1i,1.86+2.32*1i,1.9+.15*1i,2.46+.15*1i,2.46+.34*1i,2.59+.58*1i,2.72+1.51*1i,2.91+1.51*1i,3.04+.62*1i,4+.61*1i,4-.16*1i,-4-.16*1i,-4+.15*1i,-2.8+.15*1i,-2.8+.85*1i,-4+.87*1i,-4+1.54*1i];
        C = -[3.29+1.83*1i,3.45+1.69*1i,3.45+.46*1i,3.29+.28*1i,2.77+.27*1i,2.59+.46*1i,2.59+.69*1i,2.86+.69*1i,2.86+.61*1i,3.16+.61*1i,3.16+1.55*1i,2.86+1.55*1i,2.86+1.43*1i,2.59+1.43*1i,2.59+1.69*1i,2.77+1.83*1i,3.29+1.83*1i];
        O = -[2.19+1.83*1i,2.38+1.67*1i,2.38+.46*1i,2.19+.29*1i,1.62+.29*1i,1.45+.46*1i, 1.45+1.67*1i,1.62+1.83*1i, 2.19+ 1.83*1i];
        R = -[.92+1.83*1i, 1.22 + 1.83*1i, 1.22+.3*1i,.53+.3*1i,.36+.49*1i, .36+1.03*1i,.53+1.22*1i,.31+1.83*1i,.63+1.83*1i,.85+1.24*1i,.92+1.24*1i,.92+1.83*1i];
        N = -[.08+1.83*1i, .08+.3*1i,-.25+.3*1i,-.53+.96*1i,-.53+.3*1i,-.83+.3*1i,-.82+1.83*1i,-.53+1.83*1i,-.25+1.22*1i,-.25+1.83*1i,.08+1.83*1i];
        E =  [1.05-1.83*1i, 1.05-.3*1i,1.8-.3*1i,1.8-.62*1i,1.38-.62*1i,1.38-.91*1i,1.74-.91*1i,1.74-1.26*1i,1.38-1.26*1i,1.38-1.57*1i,1.81-1.57*1i,1.81-1.83*1i,1.05-1.83*1i];
        L1 = [2.04-1.83*1i,2.04-.3*1i,2.34-.3*1i,2.34-1.59*1i,2.78-1.59*1i,2.78-1.83*1i,2.04-1.83*1i];
        L2 = L1 + .93;
        
        all = {Houzes,C,O,R, N,E,L1,L2};
        all_l = length(Houzes)+length(C)+length(O)+length(R)+length(N)+length(E)+length(L1)+length(L2) - length(all)+1;
        t_l = 2*pi/(all_l-1);
        t = chebfun('t', [0,t_l], 'splitting', 'on');
        l = @(t,a,b) a + t/t_l*(b-a);
        gam1 = l(t,Houzes(1),Houzes(2));
        
        for i1 = 1:length(all)
            piec = all{i1};
            if i1 == 1
                st = 2;
            else
                st = 1;
            end
            for i2 = st:(length(piec)-1)
                gam1 = join(gam1,l(t,piec(i2),piec(i2+1)));
            end
        end
        tf = chebfun('t', [0,2*pi], 'splitting', 'on');
        gam = gam1(2*pi-tf);
    case 'HouseHolder_sky'
        Houzes = [-4 + 1.54*1i,-2.8+1.54*1i,-2.8+2.03*1i,-2.35+2.03*1i,-2.35+.22*1i,-2.26+.22*1i,-2.26+2.05*1i,-1.46+2.05*1i,-1.46+.15*1i,-1.19+.15*1i,-1.19+.53*1i,-.82+.53*1i,-.49+.89*1i,-.49+.97*1i,-.31+.97*1i,-.31+.89*1i, -.01+.53*1i,1.26+.53*1i, 1.31 + 2.32*1i,1.58+2.95*1i,1.86+2.32*1i,1.9+.15*1i,2.46+.15*1i,2.46+.34*1i,2.59+.58*1i,2.72+1.51*1i,2.91+1.51*1i,3.04+.62*1i,4+.61*1i,4-.16*1i,-4-.16*1i,-4+.15*1i,-2.8+.15*1i,-2.8+.85*1i,-4+.87*1i,-4+1.54*1i];        
        all = {Houzes};
        all_l = length(Houzes) - length(all)+1;
        t_l = 2*pi/(all_l-1);
        t = chebfun('t', [0,t_l], 'splitting', 'on');
        l = @(t,a,b) a + t/t_l*(b-a);
        gam1 = l(t,Houzes(1),Houzes(2));
        
        for i1 = 1:length(all)
            piec = all{i1};
            if i1 == 1
                st = 2;
            else
                st = 1;
            end
            for i2 = st:(length(piec)-1)
                gam1 = join(gam1,l(t,piec(i2),piec(i2+1)));
            end
        end
        tf = chebfun('t', [0,2*pi], 'splitting', 'on');
        gam = gam1(2*pi-tf);

%% 3D        
    case 'Torus'

        r = bnd_type{2};
        R = bnd_type{3};
        s = Tori(r,R);
        return

    case 'Wobbly torus'

        m = bnd_type{2}; % petal number of generating curve
        n = bnd_type{3}; % twist number along toroidal direction
        a = bnd_type{4};
        s = wobblytorus2(m,n,a,1);
        return

    case 'C Torus'

        a = bnd_type{2};
        b = bnd_type{3};
        c = bnd_type{4};
        d = bnd_type{5};
        ee = bnd_type{6};
        s = Ctor(a,b,c,d,ee);
        return
        
%% Non-smooth      
    case 'Teardrop'

        m = bnd_type{2};
        gam = chebfun(@(t) 2*((sin(2*pi-t) .* (sin((2*pi-t)/2)).^m)+cos(2*pi-t)*1i), [0, 2*pi]);

    case 'Keyhole'
        r = bnd_type{2}; R = bnd_type{3}; e = bnd_type{4};
        theta = bnd_type{5};
        t = chebfun('t',[0 2*pi/4],'splitting','on');                   % parameter
        c = [-R+e*1i -r+e*1i -r-e*1i -R-e*1i];
        gam = join( exp(1i*(pi+theta)) *(c(1) + (t/(2*pi)*4)*(c(2)-c(1))), ...       % top of the keyhole
            exp(1i*(pi+theta)) *(c(2)*c(3).^(t/(2*pi)*4) ./ c(2).^(t/(2*pi)*4)), ...    % inner circle
            exp(1i*(pi+theta)) *(c(3) + (t/(2*pi)*4)*(c(4)-c(3))), ...       % bottom of the keyhole
            exp(1i*(pi+theta)) *(c(4)*c(1).^(t/(2*pi)*4) ./ c(4).^(t/(2*pi)*4)));       % outer circle

    case 'Radiator'

        % Define parameters
        r = bnd_type{2};
        R = bnd_type{3};
        e = bnd_type{4};

        outr = bnd_type{5};
        deltatheta = bnd_type{6};

        n = bnd_type{7}-1;  % number of smaller circles

        Btheta = bnd_type{8};

        % Define inner smaller circle centers
        thetastilde = 2 * (0:n) * pi / (n+1);               % Angles for center points
        thetas = [thetastilde(thetastilde >= pi),thetastilde(thetastilde < pi)];

        centers = (R + r) / 2 * exp(1i * thetas);      % Center points for each smaller circle

        % Remove center point at -(R+r)/2
        thetas(real(centers) == -(R + r) / 2) = [];
        centers(real(centers) == -(R + r) / 2) = [];

        n = length(centers);

        % Define main keyhole points
        c = [-R - e*1i, -r - e*1i];

        for i1 = 1 : n
            theta = thetas(i1);
            cent = centers(i1);
            c = [c,r*exp((theta-deltatheta)*1i),cent + outr * exp((theta+r/outr*deltatheta+pi)*1i),cent + outr * exp((theta-r/outr*deltatheta+pi)*1i),r*exp((theta+deltatheta)*1i)];
        end

        c = [c, -r + e*1i, -R + e*1i];

        NN = length(c);

        a = 2*pi/NN;

        % Define parameter for curve
        t = chebfun('t', [0, a], 'splitting', 'on');

        gam1 = exp(1i*(pi+Btheta)) *(c(1) + (t/a)*(c(2)-c(1)));

        for i1 = 2:4:(length(c)-4)
            gam1 = join(gam1, exp(1i*(pi+Btheta)) *(c(i1)*c(i1+1).^(t/a) ./ c(i1).^(t/a)));
            gam1 = join(gam1, exp(1i*(pi+Btheta)) *(c(i1+1) + (t/a)*(c(i1+2)-c(i1+1))));

            cent = centers(ceil(i1/4));

            angle_start = angle(c(i1+2) - cent);
            angle_end = angle(c(i1+3) - cent);

            if abs(angle_end-angle_start)<pi
                angle_end = angle_end + 2*pi;
            end
            angles = [angle_start,angle_end];

            gam1 = join(gam1, exp(1i*(pi+Btheta)) *(cent + outr * exp(1i * ((1 - t/a) * angles(1) + (t/a) * angles(2)))));

            gam1 = join(gam1, exp(1i*(pi+Btheta)) *(c(i1+3) + (t/a)*(c(i1+4)-c(i1+3))));

        end

        gam1 = join(gam1, exp(1i*(pi+Btheta)) *(c(end-2)*c(end-1).^(t/a) ./ c(end-2).^(t/a)));
        gam1 = join(gam1, exp(1i*(pi+Btheta)) *(c(end-1) + (t/a)*(c(end)-c(end-1))));
        gam1 = join(gam1, exp(1i*(pi+Btheta)) *(c(end)*c(1).^(t/a) ./ c(end).^(t/a)));
        gam = chebfun(@(t)gam1(2*pi-t), [0, 2*pi], 'splitting', 'on');
    case 'Spikes'
        a = bnd_type{2};
        ts = linspace(0,2*pi,a+2);
        t = chebfun('t',[0,ts(2)-ts(1)],'splitting','on');
        spik = @(t) 2*(exp(1i*t) + 1/(a*exp(a*1i*t)));
        gam = spik(t);
        for i1 = 2:(length(ts)-1)
            gam = join(gam,spik(t+ts(i1)));
        end
    case 'Von Koch Snowflake'
        N = bnd_type{2};
        gam = VK_snowflake(N);
    case 'Perforated Wall'
        h_space = bnd_type{2};
        w_space = bnd_type{3};
        theta = bnd_type{4};
        Nholes = bnd_type{5};

        v_vert = linspace(-h_space,h_space,2*(Nholes+1));
        tot_vert = [1i*v_vert-w_space;1i*v_vert+w_space];
        tot_vert = exp(1i*theta)*tot_vert(:);
    
        t = chebfun('t',[0,2*pi/(4*(Nholes+1))]);
        a = 1/(2*pi/(4*(Nholes+1)));
        gam = tot_vert(2) * t * a+ tot_vert(1)*(1-t* a);
        gam = join(gam,tot_vert(4) * t * a + tot_vert(2)*(1-t * a));
        gam = join(gam,tot_vert(3) * t * a + tot_vert(4)*(1-t * a));
        gam = join(gam,tot_vert(1) * t * a + tot_vert(3)*(1-t * a));
    
        for i1 = 5:4:length(tot_vert)
            gam = join(gam,tot_vert(i1+1) * t * a+ tot_vert(i1)*(1-t* a));
            gam = join(gam,tot_vert(i1+3) * t * a + tot_vert(i1+1)*(1-t * a));
            gam = join(gam,tot_vert(i1+2) * t * a + tot_vert(i1+3)*(1-t * a));
            gam = join(gam,tot_vert(i1) * t * a + tot_vert(i1+2)*(1-t * a));
        end

%% Costum        
    otherwise

        load(bnd_type{1});
        if ~exist('gam')
            gam = bd;
        end

end

dgam = diff(gam);
dgam2 = diff(dgam);

s.Z = gam;
s.Zp = dgam;
s.Zpp = dgam2;
end

function s = Tori(r,R)

xf = @(u,v) (r*cos(u)+R).*cos(v);
yf = @(u,v) (r*cos(u)+R).*sin(v);
zf = @(u,v) r*sin(u);

xu =  @(u,v) -r*sin(u).*cos(v);
xv =  @(u,v) -(r*cos(u)+R).*sin(v);
xuu = @(u,v) -r*cos(u).*cos(v);
xuv = @(u,v)  r*sin(u).*sin(v);
xvv = @(u,v) -(r*cos(u)+R).*cos(v);

yu =  @(u,v) -r*sin(u).*sin(v);
yv =  @(u,v)  (r*cos(u)+R).*cos(v);
yuu = @(u,v) -r*cos(u).*sin(v);
yuv = @(u,v) -r*sin(u).*cos(v);
yvv = @(u,v) -(r*cos(u)+R).*sin(v);

zu =  @(u,v)  r*cos(u) + (v.^0-1);
zv =  @(u,v)  (u+v).^0-1;
zuu = @(u,v) -r*sin(u) + (v.^0-1);
zuv = @(u,v)  (u+v).^0-1;
zvv = @(u,v)  (u+v).^0-1;


s.Z  = @(u,v) [xf(u,v); yf(u,v); zf(u,v)];
s.Zu = @(u,v) [xu(u,v);yu(u,v);zu(u,v)];
s.Zv = @(u,v) [xv(u,v);yv(u,v);zv(u,v)];
s.Zuu= @(u,v) [xuu(u,v);yuu(u,v);zuu(u,v)];
s.Zuv= @(u,v) [xuv(u,v);yuv(u,v);zuv(u,v)];
s.Zvv= @(u,v) [xvv(u,v);yvv(u,v);zvv(u,v)];
s.topo = 't';
s.type = 'periodic'; s.L = 0;

end

function s = Ctor(a,b,c,d,ee)
    % Original curve components
    r = @(u) 3 + c * tanh(a * cos(u));
    Dr =  @(u) -a * c * sech(a * cos(u)).^2 .* sin(u);
    DDr =  @(u) -a * c * sech(a * cos(u)).^2 .* (2 * a * sin(u).^2 * tanh(a * cos(u)) + cos(u));

    th =  @(u) b * sin(u);
    Dth =  @(u) b * cos(u);
    DDth =  @(u) -b * sin(u);
    
    gam =  @(u) d * exp(1i * th(u)) .* r(u) / 2.;
    Dgam =  @(u) d * exp(1i * th(u)) .* (1i * Dth(u) .* r(u) + Dr(u)) / 2;
    DDgam =  @(u) d * exp(1i * th(u)) .* (1i * Dth(u) .* r(u) + Dr(u)) .* 1i .* Dth(u) - exp(1i * th(u)) .* (1i * DDth(u) .* r(u) + 1i * Dth(u) .* Dr(u) + DDr(u))/2;

    
    % Extract x and y components from gam
    x = @(u) real(gam(u)) + ee;
    z = @(u) imag(gam(u));
    Dx = @(u) real(Dgam(u));
    Dz = @(u) imag(Dgam(u));
    DDx = @(u) real(DDgam(u));
    DDz = @(u) imag(DDgam(u));

    xf = @(u,v) cos(v) .* x(u);
    yf = @(u,v) sin(v) .* x(u);
    zf = @(u,v) z(u) + v.^0 - 1;

    % [u, v] = meshgrid(linspace(0,2*pi,200),linspace(0,2*pi,200));
    % u = u(:); v = v(:);
    % plot3(xf(u,v),yf(u,v),zf(u,v),'.k');

    xu = @(u,v) cos(v) .* Dx(u);
    xv = @(u,v) -sin(v) .* x(u);
    xuu = @(u,v) cos(v) .* DDx(u);
    xuv = @(u,v) -sin(v) .* Dx(u);
    xvv = @(u,v) -cos(v) .* x(u);

    yu = @(u,v) sin(v) .* Dx(u);
    yv = @(u,v) cos(v) .* x(u);
    yuu = @(u,v) sin(v) .* DDx(u);
    yuv = @(u,v) cos(v) .* Dx(u);
    yvv = @(u,v) -sin(v) .* x(u);

    zu = @(u,v) Dz(u) + v.^0 - 1;
    zv = @(u,v) (v+u).^0 - 1;
    zuu = @(u,v) DDz(u) + v.^0 - 1;
    zuv = @(u,v) (v+u).^0 - 1;
    zvv = @(u,v) (v+u).^0 - 1;

    s.Z  = @(u,v) [xf(u,v); yf(u,v); zf(u,v)];
    s.Zu = @(u,v) [xu(u,v);yu(u,v);zu(u,v)];
    s.Zv = @(u,v) [xv(u,v);yv(u,v);zv(u,v)];
    s.Zuu= @(u,v) [xuu(u,v);yuu(u,v);zuu(u,v)];
    s.Zuv= @(u,v) [xuv(u,v);yuv(u,v);zuv(u,v)];
    s.Zvv= @(u,v) [xvv(u,v);yvv(u,v);zvv(u,v)];
    s.topo = 't';
    s.type = 'periodic'; s.L = 0;

end

function gam = VK_snowflake(N)
p2 = 1; 
p1 = 1/2 + 1i * sqrt(3)/2;            
p0 = 0;                    

% Get final points from recursion
final_points = koch_snowflake([p0, p1, p2, p0], N);

final_points = 6.5*(final_points - mean(final_points));
num_pts = length(final_points);

t = chebfun('t',[0,2*pi/(num_pts-1)]);
gam = final_points(2)*(t)/(2*pi/(num_pts-1)) + final_points(1)*(1-t/(2*pi/(num_pts-1)));
for i1 = 2:(num_pts-1)
    gam = join(gam, final_points(i1+1)*(t)/(2*pi/(num_pts-1)) + final_points(i1)*(1-t/(2*pi/(num_pts-1))));
end

end
% Function to generate final points of Koch curve
function final_points = koch_snowflake(p, n)
    if n == 0
        final_points = p; % Return points at base case
    else
        new_p = [];
        for k = 1:length(p)-1
            % Divide segment into three parts
            z1 = p(k);
            z2 = p(k+1);
            v = (z2 - z1) / 3;
            
            % Compute new points
            a = z1 + v;
            b = z1 + 2*v;
            c = a + v * exp(1i * pi/3); % Rotate by 60 degrees

            % Store new sequence
            new_p = [new_p, z1, a, c, b]; 
        end
        new_p = [new_p, p(end)]; % Add last point
        
        % Recursively refine
        final_points = koch_snowflake(new_p, n-1);
    end
end



