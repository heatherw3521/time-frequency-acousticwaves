function s = ns_setupquad(s, Ninp, plot_me)
% NS_SETUPQUAD  
% This code is taken from Alex Barnett 10/8/14 and adapted to attempt a
% correction scheme for panel based quadrature. It gives visually correct
% results, but mathematically is incorrect. It must be replaced by
% something like Kolm-Rohklin quadrature.
% Inputs:
%  s : struct containing either the field s.Z, a mapping from [0,2pi) to
%      the complex plane, which must be able to accept vectorized inputs,
%      or the field s.x containing column vector of the nodes.
%  N : number of nodes (not needed if s.x is present; however, if s.Z is also
%      present, N overrides the number of nodes in s.x).
% Outputs:
%  s : same struct with added column vector fields (nodes, weights, velocities,
%      curvatures, etc) needed for quadrature and Nystrom methods. Namely,
%      s.x nodes in complex plane, Z(s_j)
%      s.xp velocities Z'(s_j)
%      s.xp accelerations Z''(s_j)
%      s.t parameter values s_j for trapezoid rule, 2pi.(j-1)/N, j=1,...,N
%      s.nx outward unit normals
%      s.tang forward unit tangential vectors
%      s.sp speeds |Z'(s_j)|
%      s.w "speed weights" (2pi/N)*s.sp
%      s.cw velocity weights (ie complex speed)
%      s.cur curvatures kappa(s_j) (inverse bending radius)
%
% Example usage:
%
% s.Z = @(s) (1+0.3*cos(5*s).*exp(1i*s);                   % starfish param
% s = setupquad(s,100);
% figure; plot(s.x,'k.'); hold on; plot([s.x, s.x+0.2*s.nx].', 'b-'); axis equal
%
% Now check that normals from spectral differentiation are accurate:
%
% s.Zp = @(s) -1.5*sin(5*s).*exp(1i*s) + 1i*s.Z(s);        % Z' formula
% t = setupquad(s,100); norm(t.nx-s.nx)                 % should be small

% (c) Alex Barnett 10/8/14, name changed to avoid conflict w/ mpspack 6/12/16.
% 0-indexed to match interp, 6/29/16


domdom = s.Z.domain;

if nargin>1           % use N from args
    %non-smooth point
    if length(domdom) == 2
        if abs(s.Zp(0)) < 1e-13
            tt_point = [0,2*pi];
        else

            curvi = abs(real(s.Zp) * imag(s.Zpp) - imag(s.Zp)*real(s.Zpp) ) / (real(s.Zp)^2+imag(s.Zp)^2)^(3/2);
            Dcurvi = diff(curvi);
            DDcurvi = diff(Dcurvi);
            
            tt_point = roots(Dcurvi);
            tt_point = tt_point(DDcurvi(tt_point)<-.1);
            if length(tt_point)>4
                tt_point = tt_point(abs(diff(tt_point))>.1);
                tt_point = tt_point(abs(tt_point)>.1);
            end
            if length(tt_point) == 2
                tt_point = [tt_point.',tt_point(1)+2*pi];
            else
                tt_point = [tt_point(end)-2*pi,tt_point.',tt_point(1)+2*pi];
            end
        end
    elseif abs(s.Z(s.Z.domain(2))) < 1e-13

        NN = length(s.Z.domain);
        all_tt_point = [];
        for i3 = 1:(NN-1)
            Z = chebfun(s.Z, [domdom(i3)+1e-15,domdom(i3+1)-1e-15]);
            Zp = diff(Z);
            Zpp = diff(Zp);
            curvi = abs(real(Zp) * imag(Zpp) - imag(Zp) * real(Zpp) ) / (real(Zp)^2+imag(Zp)^2)^(3/2);
            Dcurvi = diff(curvi);
            DDcurvi = diff(Dcurvi);
            
            tt_point = roots(Dcurvi);
            tt_point = tt_point(DDcurvi(tt_point)<-1);
            if length(tt_point)>4
                iindi = find(diff(tt_point)>.15);
                for i4 = 1:length(iindi)
                    if iindi(i4) > length(tt_point)/2
                        iindi(i4) = iindi(i4) + 1;
                    end
                end
                tt_point = tt_point(iindi);
                tt_point = tt_point(abs(tt_point)>.1);
            end
            tt_point = tt_point(abs(Dcurvi(tt_point))<1e-5);
            all_tt_point = [all_tt_point; tt_point];
        end
        tt_point = all_tt_point;
        if length(tt_point) == 2
            tt_point = [tt_point.',tt_point(1)+2*pi];
        elseif isempty(tt_point)
            tt_point = s.Z.domain;
        else
            tt_point = [tt_point(end)-2*pi,tt_point.',tt_point(1)+2*pi];
        end
    else
        tt_point = s.Z.domain;
    
    end

    %panel points
    NN = Ninp{2};

    tt_point = tt_point(tt_point > 1e-10 & tt_point < 2*pi - 1e-10);
    tt_point = [0, tt_point, 2*pi];

    tt_point = panel_halver(s, tt_point, NN);
    

    n = Ninp{1};
    t = [];
    w = [];
    n_adp = adp_points_nr(n,s, tt_point);
    if any(n_adp==0)
       n_adp = n_adp+1; 
    end
    for i1 = 1:(length(tt_point)-1)
        [tt,ww] = lgwt(n_adp(i1), tt_point(i1), tt_point(i1+1));
        %[tt,ww] = lgwt_points_on_curve(n, tt_point(i1), tt_point(i1+1), s);
        
        t = [t;tt];
        w = [w;ww];
    end
    [t, ia] = unique(mod(t,2*pi), 'stable');
    w = w(ia);
    [sorted_t, sort_idx] = sort(t);
    sorted_w = w(sort_idx);
    diffs = [true; diff(sorted_t) > 1e-12];
    t = sorted_t(diffs);
    w = sorted_w(diffs);
    s.t = t;
    if isfield(s,'Z'), s.x = s.Z(s.t); end    % use formula
    if plot_me
        figure(23)
        plot_pannel(s,mod(tt_point,2*pi),s.t); title('Nystrom discretization','FontSize',20)
    end
    %if Ntot~=length(s.x), error('N differs from length of s.x; that sucks!'); end
else
    error('Need to provide at least s.Z and N, or s.x. Neither found!');
end
if isfield(s,'Zp'), s.xp = s.Zp(s.t); else, s.xp = perispecdiff(s.x); end
if isfield(s,'Zpp'), s.xpp = s.Zpp(s.t); else, s.xpp = perispecdiff(s.xp); end
% Now local stuff that derives from x, xp, xpp at each node...
s.sp = abs(s.xp);
s.tang = s.xp./s.sp;
s.nx = -1i*s.tang;
s.cur = -real(conj(s.xpp).*s.nx)./s.sp.^2;  % recall real(conj(a)*b) = "a dot b"
s.w = w.*s.sp;
s.cw = w.*s.xp;  % complex weights (incl complex speed)

function g = perispecdiff(f)
% PERISPECDIFF - use FFT to take periodic spectral differentiation of vector
%
% g = perispecdiff(f) returns g the derivative of the spectral interpolant
%  of f, which is assumed to be the values of a smooth 2pi-periodic function
%  at the N gridpoints 2.pi.j/N, for j=1,..,N (or any translation of such
%  points). Can be row or col vec, and output is same shape.
%
% Without arguments, does a self-test.

% Barnett 2/18/14
N = numel(f);
if mod(N,2)==0   % even
    g = ifft(fft(f(:)).*[0 1i*(1:N/2-1) 0 1i*(-N/2+1:-1)].');
else
    g = ifft(fft(f(:)).*[0 1i*(1:(N-1)/2) 1i*((1-N)/2:-1)].');
end
g = reshape(g,size(f));





function [tt_point] = panel_halver(s, tt_point, NN)

if NN == 0
    return
end
for i1 = length(tt_point):-1:2

    gam = chebfun(@(t) s.Z(mod(t,2*pi)), [tt_point(i1-1),tt_point(i1)], 'splitting', 'on');

    len = cumsum(abs(diff(gam)));
    
    L = arcLength(gam);

    tt_temp_point = roots(len - L/2);


    for i2 = 2:NN
        gam = chebfun(@(t) s.Z(mod(t,2*pi)), [tt_point(i1-1),tt_temp_point(1)], 'splitting', 'on');
        L = arcLength(gam);
        len = cumsum(abs(diff(gam)));
        rts_L2 = (roots(len - L/2));
        tt_temp_point = [rts_L2(1),tt_temp_point];

        gam = chebfun(@(t) s.Z(mod(t,2*pi)), [tt_temp_point(end),tt_point(i1)], 'splitting', 'on');
        L = arcLength(gam);
        len = cumsum(abs(diff(gam)));
        tt_temp_point = [tt_temp_point, roots(len - L/2)];
    end
    tt_point = [tt_point(1:i1-1),tt_temp_point,tt_point(i1:end)];
end


function [t,w] = lgwt_points_on_curve(n, a, b, s)

len = cumsum(abs(diff(s.Z)));

[st, w] = lgwt(n,len(a),len(b));


t = st*0;
for i1 = 1:length(st)
    rts = roots(len-st(i1));
    if length(rts) == 2
        t(i1) = rts(1);
    else
        t(i1) = rts;
    end

end



function [] = plot_pannel(s,pannel_t,t)

plot(s.Z, 'k', 'linewidth',1); hold on;
mydisc_fill(real(s.Z(linspace(0,2*pi,8001))), imag(s.Z(linspace(0,2*pi,8001))), [.7 .7 .7]);

normal  = (-1i * s.Zp(pannel_t))./abs(s.Zp(pannel_t));
normal_length = 0.1;

points = s.Z(pannel_t);

start_pt = points - normal * normal_length;
end_pt = points + normal * normal_length;
for i1 = 1:length(start_pt)
    plot(real([start_pt(i1), end_pt(i1)]),imag([start_pt(i1), end_pt(i1)]), 'r', 'LineWidth', 2);
end
scatter(real(s.Z(t)),imag(s.Z(t)),'*k')
axis square
set(gca,'visible','off')

function h = mydisc_fill(x, y, clr)
n = length(x);
dz = sqrt((x(2:n) - x(1:(n-1))).^2+(y(2:n) - y(1:(n-1))).^2);

ex_I = find(dz>=.5);
plt_points = [0;ex_I(:);n];

hold on
for i1 = 1:(length(plt_points)-1)
    indi = [(plt_points(i1)+1):plt_points(i1+1),(plt_points(i1)+1)];
    h = fill(x(indi), y(indi),clr);
    %h.FaceAlpha = .75; %some transparency
end

function n_adp = adp_points_nr(n, s, t)

n_adp = zeros(length(t)-1,1);
for i1 = 1:(length(t)-1)
    t1 = t(i1); t2 = t(i1+1);

    L = arcLength(s.Z);
    l = arcLength(chebfun(s.Z,[t1 + 1e-15,t2 - 1e-15])); 
    n_adp(i1) = floor(n * l / L);
end
allt_i = 1:length(t)-1;
while sum(n_adp) ~= n
    min_in = allt_i(n_adp == min(n_adp));
    n_adp(min_in(1)) = n_adp(min_in(1)) + 1; % Can be made a lotm ore efficiently
end

