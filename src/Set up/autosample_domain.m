function [zz, xx, yy, keep] = autosample_domain(n, rnge, bnd_type, plotme)
xrng = rnge{1}; yrng = rnge{2};
x = linspace(xrng(1), xrng(2), n).';
y = linspace(yrng(1), yrng(2), n);
[xx,yy] = meshgrid(x.', y);
zz = xx + 1i*yy;
sz = size(zz);
zz = zz(:);
s = setup_boundary(bnd_type);
gam = s.Z;
%now keep only the points exterior to the object
if isa(gam, 'chebfun') && length(s.Z.domain) == 2
    keep = ~incontour(gam, {zz});
elseif length(rnge) == 2
    switch bnd_type{1}
        case 'Circles'
            c = bnd_type{2}; % centers
            R = bnd_type{3}; %radii
            keep = abs(zz-c(1))<=R(1);
            for i1 = 2:length(R)
                keep = keep | abs(zz-c(i1))<=R(i1);
            end
            keep = ~keep;
        case 'Whisper Gallery'
            a = bnd_type{2}; %
            b = bnd_type{3};
            fattener = bnd_type{4};

            keep = ~ ((real(zz)/a).^2+(imag(zz)/b).^2 >= 1 & (real(zz)/(a*fattener)).^2+(imag(zz)/(b/fattener)).^2 <= 1);
        case 'Crescent Gallery'
            r = 5; 
            th = exp(1i*pi/2);
            p = 1i; 
            a1 = bnd_type{2}; % ?
            a2 = bnd_type{3}; % ?
            dist = bnd_type{4}; % Distance between crescents

            gam1 = r*th*(chebfun(@(t) p+conj(exp(1i*t)- a1/(exp(1i*t)+a2))+dist/2, [0, 2*pi]));
            gam2 = r*th*chebfun(@(t) -conj(exp(1i*t)- a1/(exp(1i*t)+a2))-dist/2, [0, 2*pi]);
            keep = ~(incontour(gam1, {zz(:)}) | incontour(gam2, {zz(:)}));
        case 'Perforated Wall'
        h_space = bnd_type{2};
        w_space = bnd_type{3};
        theta = bnd_type{4};
        Nholes = bnd_type{5};

        v_vert = linspace(-h_space,h_space,2*(Nholes+1));
        
        keep = zz == 0;

        Mz = max(abs(zz));
        for i1 = 1:2:length(v_vert)
            keep = keep | (imag(exp(-1i*theta)*zz) > v_vert(i1) - 1e-3 / Mz  & v_vert(i1+1) + 1e-3 / Mz > imag(exp(-1i*theta)*zz) & abs(real(exp(-1i*theta)*zz)) < w_space + 1e-3);
        end
        
        keep = ~keep;

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
        


        nkeep = zeros(size(zz));
        for i1 = 1:length(all)
            
            piec = all{i1};

            t_l = 2*pi/(length(piec)-1);
            t = chebfun('t', [0,t_l], 'splitting', 'on');
            l = @(t,a,b) a + t/t_l*(b-a);

            gam1 = l(t,piec(1),piec(2));

            for i2 = 2:(length(piec)-1)
                gam1 = join(gam1,l(t,piec(i2),piec(i2+1)));
            end
            tf = chebfun('t',[eps,2*pi-eps],'splitting','on');
            gam = gam1(2*pi-tf);
            nkeep = nkeep | incontour(gam, {zz});

        end
        keep = ~nkeep;
        otherwise
            if isa(gam, 'chebfun')
                keep = ~incontour(gam, {zz});
            end
    end
    
else
    zrng = rnge{3};

    [XX,YY] = meshgrid(linspace(xrng(1),xrng(2),n),linspace(yrng(1),yrng(2),n));
    xx = repmat(reshape(XX,1,[]),1,n);
    yy = repmat(reshape(YY,1,[]),1,n);

    ZZ = linspace(zrng(1),zrng(2),n);
    zz = []; % This can be made better
    for i = 1:n
        repeated_elements = repmat(ZZ(i), 1, n^2);
        zz = [zz, repeated_elements];
    end
    switch bnd_type{1}
        case 'Torus'

            r = bnd_type{2};
            R = bnd_type{3};

            keep = ((R-sqrt(xx.^2+yy.^2)).^2+zz.^2-r^2) > 0;

        otherwise
            keep = ~incontour(gam, {xx,yy,zz});
            if strcmp(bnd_type{1},'C Torus')
                zuni = unique(zz);
                nn = sum(zuni(1) == zz);
                ntot = length(zz)/nn;
                for i1 = 1:nn
                    if ~(sum(keep(i1+(0:(ntot-1))*nn))== length(keep(i1+(0:(ntot-1))*nn)))
                        keep(i1+(0:(ntot-1))*nn) = ~keep(i1+(0:(ntot-1))*nn);
                    end
                end
            end
    end
end

if plotme
    figure(3);
    if bnd_type{1} == "Torus" || bnd_type{1} == "Wobbly torus" || bnd_type{1} == "C Torus"
        plot3(xx(~keep),yy(~keep),zz(~keep), 'x'); title("Points inside domain",'FontSize',20)
        xlabel("x",'FontSize',20);ylabel("y",'FontSize',20);zlabel("z",'FontSize',20);
    else
        plot(zz(keep), 'x'); title("Points outside domain",'FontSize',20)
        xlabel("x",'FontSize',20);ylabel("y",'FontSize',20);
    end
end
end


function out = incontour(gamma, xs)
% determine whether x is inside gamma (y = 1), or outside(y=0).
% here gamma is a closed contour.

if isscalar(xs)
    x = xs{1};
    t = linspace(0, 2*pi, 501);
    t = t(1:end-1).';
    %N = (length(t)+1)/2;
    %sz = size(x);
    % to avoid sampling too close to the boundary, we "push out" the 
    % boundary in the normal direction. There will be issues here
    % if the boundary has corners or cusps, or potentially related to
    % nonconvexity. For now, it should work.

    n = normal(gamma);
    z = gamma(t); nn = n(t); nn = nn./sqrt(nn(:,1).^2 + nn(:,2).^2);
    nn = nn(:,1) + 1i*nn(:,2);
    zz = z + .05*nn;
    out = inpolygonc(x,zz);
else
    xx = xs{1}; yy = xs{2}; zz = xs{3};
    NN = 100;
    theta = linspace(0, 2*pi, NN);
    phi = linspace(0, 2*pi, NN);

    [Theta, Phi] = meshgrid(theta, phi);

    X = Theta*0; Y = Theta*0; Z = Theta*0;

    for i1 = 1:NN
        for i2 = 1:NN
            G = gamma(Theta(i1,i2),Phi(i1,i2));
            X(i1,i2) = G(1); Y(i1,i2) = G(2); Z(i1,i2) = G(3);
        end
    end

    fv = surf2patch(X, Y, Z, 'triangles');
    vertices = fv.vertices;
    faces = fv.faces;

    out = inpolyhedron(faces, vertices, [xx',yy',zz']);
end
%dgamma = diff(gamma);
%y = gamma(t);
%out = sum(pi/N*(1./(y-x.').*dgamma(t)));
%out = out(:);
%out(abs(out)>1e-5) = 1; %in contour
%out(abs(out)<=1e-5) = 0; %outsided contour
%out = logical(out);
end

function T = inpolygonc(z,w) %complex version of inpolygon
[T, ON] = inpolygon(real(z),imag(z),real(w),imag(w));
%correct so that points on the edge are not included:
T(ON) = false;
end