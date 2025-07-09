function [uinc_res, Uinc_res, Uinc_slow_res, TT, W, windows] = setup_incident(wave_type, comp_type_freq, dimension, s)

switch wave_type{1}
    case 'Tsunami'

        %Initializing parameters:
        w0 = wave_type{2};
        c = wave_type{3};
        sigma = wave_type{4};
        rs = wave_type{5};
        
        T = wave_type{6};
        TT = 2 * T; %TT = 60;

        %define incident field in frequency space:
        uinc_part = @(t,A) 3/(sqrt(2*pi)*sigma) * exp(-1i*w0*(t-A)) .* exp(-(t-A).^2/(2*sigma^2));
        Uinc_part = @(w,A) 3*exp(-(w-w0).^2*sigma^2/2) .* exp(1i*A.*w);

        if dimension == 2
            Uinc = @(z, w) 0;
            uinc = @(z, t) 0;
            for i1 = 1:length(rs)
            
                r = rs(i1);
                r = r/abs(r); %incident field direction of prop.
                Uinc = @(z, w) Uinc(z,w) + Uinc_part(w,real(r).*real(z)/c + imag(r).*imag(z)/c+T);
                uinc = @(z, t) uinc(z,t) + uinc_part(t,real(r).*real(z)/c + imag(r).*imag(z)/c+T);
            end
        elseif dimension == 3
            r = rs/sqrt(sum(rs.^2)); %incident field direction of prop.
            Uinc = @(z, w) Uinc_part(w,r(1).*z(1,:)/c + r(2).*z(2,:)/c + r(3).*z(3,:)/c+T);
            uinc = @(z, t) uinc_part(t,r(1).*z(1,:)/c + r(2).*z(2,:)/c + r(3).*z(3,:)/c+T);
        end
        
    otherwise
        error("Unrecognized incident wave type")
end

% Frequency band
if dimension == 2
    scat_points = s.Z(linspace(0,2*pi,500));
else
    [thet,phi] = meshgrid(linspace(0,2*pi,50),linspace(0,2*pi,50));
    scat_points = s.Z(thet(:).',phi(:).');
end
TOL = 1e-13;
W = freq_band_finder(Uinc, scat_points, TOL);
% Windowing
[uinc_res, Uinc_res, Uinc_slow_res, windows] = incident_comp_form(comp_type_freq,uinc,Uinc,dimension,TT);

end

function [uinc, Uinc, Uinc_slow, windows] = incident_comp_form(comp_freq,uinc_full,Uinc_full,dimension,T)

comp_freq_type = comp_freq{1};
windows = [0];
Uinc_slow = {};
switch comp_freq_type
    case 'known'
        uinc = {uinc_full};
        Uinc = {Uinc_full};
    case 'basic'
        N_ft = comp_freq{2};
        uinc = {uinc_full};
        Uinc = {@(z,w) Comp_freq(uinc_full,T,z,w,N_ft)};
    case 'window all'
        %W = [-10, 10];
        N_ft = comp_freq{2};
        N_win = comp_freq{3};
        windows = linspace(0,T,N_win);
        H = (windows(2)-windows(1))/3*2;

        for i1 = 1:length(windows)
            uinc{i1} = @(z,t) Window_ft(@(z,tt) uinc_full(z,tt),z, t, windows(i1), N_ft, H, "window function");
            Uinc{i1} = @(z,w) Window_ft(@(z,t) uinc_full(z,t),z, w, windows(i1), N_ft, H);
        end

    case 'window one'

        c = 1;
        if dimension == 2
            r = 1+.5*1i;
            r = r/abs(r); %incident field direction of prop.
            Uinc = {@(z,w) exp(1i* w .* (real(r).*real(z)/c + imag(r).*imag(z)/c+T))};
        elseif dimension == 3
            r = [0,0,1];
            r = r/sqrt(sum(r.^2)); %incident field direction of prop.
            Uinc = {@(z,w) exp(1i* w .* (r(1).*z(1,:)/c + r(2).*z(2,:)/c + r(3).*z(3,:)/c+T))};
        end

        N_ft = comp_freq{2};
        N_win = comp_freq{3};
        windows = linspace(0,T,N_win);
        H = (windows(2)-windows(1))/3*2;

        for i1 = 1:length(windows)
            uinc{i1} = @(z,t)      Window_ft(@(z,tt) uinc_full(z,tt),z, t, windows(i1), N_ft, H, "window function");
            if dimension == 3
                Uinc_slow{i1} = @(z,w) Window_ft(@(z,t) uinc_full(z,t),z, w, windows(i1), N_ft, H) .* (transpose(Uinc{1}(z,w)) .^ (-1));
            else
                Uinc_slow{i1} = @(z,w) Window_ft(@(z,t) uinc_full(z,t),z, w, windows(i1), N_ft, H) .* (Uinc{1}(z,w) .^ (-1));
            end
        end

        if dimension == 2
            z1 = .23+1i;
            z2 = 2.7;
        elseif dimension == 3
            z1 = [0;0;0];
            z2 = [2.7;.1;-1];
        end

        % if sum(abs((Uinc_slow{1}(z1,rand()) - Uinc_slow{1}(z2,rand()))))>=1e-13
        %     error('The Fourier transform is not seperable. Use "window all" setting.') % Does not always work            
        % end
end
end
function uinc = uinc_bump(z,t,sigma,c,r,c0,w0,T)
uinc = zeros(length(t),size(z,2));
for i1 = 1:length(t)
     uinc(i1,:) = 1/(sqrt(2*pi)*sigma) * exp(-sum((z - (c0 + c * r * t(i1))).^2)/(2*sigma^2)) .* exp(-1i*w0*(t(i1)-(r(1).*z(1,:)/c + r(2).*z(2,:)/c + r(3).*z(3,:)/c+T))); %sum() sums rows
end
end


function u = Comp_freq(f,T,z,w,N)
I = [0,T];
u = zeros(length(z), length(w));
for i2 = 1:length(w)
    ww = w(i2);
    u(:,i2) = Clen_Curt(@(z,t) f(z,t) .* exp(1i*t*ww) ,I,N,z);
end
end


function w = window(tt,T)
w = tt * 0;
eta = @(u) exp(2*exp(-1./u)./(u-1));
for i1 = 1:length(tt)
    t = tt(i1);
    if (20< t && t< 35)
        w(i1) = 1-eta((t-20)/15);
    elseif (35<= t && t<= T-35)
        w(i1) = 1 ;
    elseif (T-35< t && t< T-20)
        w(i1) = eta((t-(T-35))/15);
    end
end
end

function W = freq_band_finder(Uinc, zz, TOL)
    %freq band: Uinc is supported on band W. It is
    % dependent on parameters above. If you change them, check the
    % support on scatterer points zz with tolerance TOL 
    % and it is chosen appropriately.
    NN = 1000;
    k = linspace(-100,100,NN);
    indi = 1:NN;
    tflist = indi == 0;
    for i1 = indi
        tflist(i1) = max(abs(Uinc(zz,k(i1)))) > TOL;
    end
    W(1) = min(k(tflist));
    W(2) = max(k(tflist));
end