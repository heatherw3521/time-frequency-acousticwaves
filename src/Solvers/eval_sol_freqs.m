function outobj = eval_sol_freqs(inobj, zz, field_type, eval_type, Freqvals)
% evaluate the HH solution U(w,zz) in frequency space
% at spatial locations zz and freqs = freqvals.
%
% Freqvals must be a subset of sol.freq (frequencies we have solved for).
% if Freqvals is empty, we evaluate over all the values in sol.freq.
% (i.e., our approx to the whole support of sol)
%
% out is a matrix with cols that vary with freq, rows
% that vary with location, location is col vector stacked.
%
% 'field_type' = 'scatter', 'total', 'incident'.
% w = frequency

T = inobj.tau;
if ~isempty(Freqvals) %frequency specified:
    w = Freqvals;
else
    w = inobj.freq;
end

%% Incident field_type case
Uinc = inobj.incF;
if strcmpi(field_type, 'incident') || strcmpi(field_type, 'total')
        if ~isa(Uinc, 'function_handle')
            error("eval_sol_freq: Incident wave cannot be discrete.")
        end
        outobj.tau = [];
        Us_inc = zeros(length(zz),length(w));
        for i2 = 1:length(w)
            Us_inc(:,i2) = Uinc(zz,w(i2));
        end
        outobj.freqvals = Us_inc;
    if strcmpi(field_type, 'incident')
        return
    end
end

s = inobj.s;

dimension = 2;
if isfield(s, 'Zu')
    dimension = 3;
end

if dimension == 2
    tg.x = zz(:); %target points
   % inc = inobj.incF;
   % Uinc = zeros(length(tg.x), length(w)); %freqs var = cols, loc var = rows.
else
    tg.x = zz; %target points
   % inc = inobj.incF;
   % Uinc = zeros(size(zz,2), length(w)); %freqs var = cols, loc var = rows.
end

%for scattered field we have several strategies:
switch eval_type
    case 'basic'
        %scattered field:
        Usc = zeros(length(tg.x), length(w));
        if dimension == 2
            J = 1:length(w);
            st = sum(abs(T)); 
            fil = find(st > 1e-13); 
            J = J(fil);
        else
            J = 1:length(w);
        end
        for j = J
        %for j = 1:length(w)
            ka = w(j);
            ieta = 1i*real(ka);
            % ieta = 1i*ka;
            tau = T(:,j);
            if dimension == 2
                if abs(ka) > 1e-16
                    A = HelmDLP(ka,tg,s)-ieta*HelmSLP(ka,tg,s);
                else
                    A = LapDLPmat(tg,s) + s.w.';
                end
            elseif dimension == 3
                if abs(ka) > 1e-16   
                    A = Helm3dDLPmat(ka,tg,s) - ieta*Helm3dSLPmat(ka,tg,s);
                else 
                    zs = [0.3,-0.9,0; 0.5,0.85,0; -0.99,-0.1,0].';    % source & strength inside torus
                    A = Lap3dDLPmat(tg,s) + s.w./vecnorm(tg.x-zs(:,1))';
                end
            end
            % A = bsxfun(@times, A, s.w(:).');
            Usc(:,j) = A*tau;
        end

    case 'fmm' % NOT IMPLEMENTED WITH LAPLACE EQ

        Usc = zeros(length(tg.x), length(w));
        if dimension == 2
            J = 1:length(w);
            st = sum(abs(T)); 
            fil = find(st > 1e-13); 
            J = J(fil);
        else
            J = 1:length(w);
        end
        for j = J
            if dimension == 2
                ka = w(j); %wave number
                tau = T(:,j);
                %params = [2*pi, ka];
                %Usc(:,j) = LOCAL_apply_HBS(nodes, tau); %matvec with tau
                Usc(:,j) = wrap_fmm_targ(zz,tau,ka, s);
            elseif dimension == 3
                ka = w(j);
                tau = T(:,j);
                Usc(:,j) = wrap_fmm3d_targ(zz, transpose(tau), ka, s);
            end
            
        end
end

if strcmpi(field_type, 'scatter')
    outobj = Usc;
elseif strcmpi(field_type, 'total')
    outobj = Usc + Us_inc;
else
    error('eval_sol_freq: field type not recognized')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uu = wrap_fmm_targ(zz,tau,ka,s)
%wrapper for calling 2Dfmm to evaluate the solution, which is given by
%A = HelmDLP(ka,tg,s)-ieta*HelmSLP(ka,tg,s); out = A*tau
%
% This function has a dependency on (1) fmm2dlib, (2) fftw packages.
iprec = 4; %1e-12 tol
nsource = length(s.x); %quad nodes = source points
source = [real(s.x).'; imag(s.x).'];
ifcharge = 1;
charge = -1i*real(ka)*s.w.*tau; % quad weights and tau and mixed layer pot. param
ifdipole = 1;
dipstr = s.w.*tau; %dipole strengths (go with double layer pot)
dipvec = [real(s.nx).'; imag(s.nx).']; %normals for the source points.
ifpot = 0;
ifgrad = 0;
ifhess = 0; % we do not want potential, etc, at sources.
ntarget = length(zz);
target = [real(zz).'; imag(zz).'];
ifpottarg = 1; % we want potential at target points.
ifgradtarg = 0;
ifhesstarg = 0;

uu =hfmm2dpart(iprec,ka,nsource,source,ifcharge,charge,ifdipole,...
    dipstr,dipvec,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,...
    ifgradtarg,ifhesstarg);

uu = uu.pottarg;
end


function uu = wrap_fmm3d_targ(zz,tau,ka,s)
%wrapper for calling fmm3d to evaluate the solution, which is given by
%A = HelmDLP(ka,tg,s)-ieta*HelmSLP(ka,tg,s); out = A*tau
%
% This function has a dependency on (1) FMM3D.
tol = 1e-12; % Precision

srcinfo.sources = s.x; % double(3,n);   
% source locations, $x_{j}$

srcinfo.charges = -1i*real(ka)*s.w .* tau; 
% charge densities, $c_{j}$ (optional, 
% default - term corresponding to charges dropped)

srcinfo.nd = 1;
%           number of charge/dipole vectors (optional, 
%           default - nd = 1)

dipstr = s.w .* tau; %?????????????? I put it with the dipoles source info
%dipole strengths (go with double layer pot)

srcinfo.dipoles(1,:,:) = s.nx .* (s.w .* tau);
% dipole orientation vectors, $v_{j}$ (optional
% default - term corresponding to dipoles dropped) 

pg = 1;
%  -  pg: integer
%        | source eval flag
%        | potential at sources evaluated if pg = 1
%        | potential and gradient at sources evaluated if pg=2

targ = zz; % target locations, $t_{i}$ 
pgt = 1; % potential at targets evaluated if pgt = 1

uu = hfmm3d(tol, ka, srcinfo, pg, targ, pgt);

uu = uu.pottarg;
end


%function uu = LOCAL_fmm_targ(Ctarg,qq,C,flag_pot,params, s)

%curvelen = params(1);
%ntot = size(C,2);
%kh = params(2);
%h = curvelen/ntot;

%nsource = size(C,2);
%ntarget = size(Ctarg,2);

%target = [Ctarg(1,:);
%          Ctarg(4,:)];

% switch flag_pot
% case {'hd'}
%   ifcharge = 0;
%   charge   = zeros(1,nsource);
% case {'ht'}
%   ifcharge = 1;
%   speed    = sqrt(C(2,:).^2+C(5,:).^2);
%   %speed = s.sp; speed = speed(:);
%   charge   = 1i*kh/ntot*(qq.').*speed;
% otherwise
%   error('Layer potential not supported.')
% end
%
% ifdipole   = 1;
% ifpot      = 0;
% ifgrad     = 0;
% ifhess     = 0;
% ifpottarg  = 1;
% ifgradtarg = 0;
% ifhesstarg = 0;
%
% iprec = 5;
%
% dipvec = zeros(2,nsource);
% %dipvec(1,:) = real(s.nx);
% %dipvec(2,:) = -imag(s.nx);
% dipvec(1,:) = C(5,:);
% dipvec(2,:) = -C(2,:);
%
% %dipvec = dipvec./vecnorm(dipvec);
%
% source      = zeros(2,nsource);
% %source(1,:) = real(s.x);
% %source(2,:) = imag(s.x);
% source(1,:) = C(1,:);
% source(2,:) = C(4,:);
%
% dipstr = qq.';
% %[U]=hfmm2dpart(iprec,zk,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,ifgrad,ifhess,ntarget,target,ifpottarg,ifgradtarg,ifhesstarg)
% [U] = hfmm2dpart(iprec,kh,nsource,source,ifcharge,charge,ifdipole,...
%                  dipstr,dipvec,ifpot,ifgrad,ifhess,ntarget,target,...
%                  ifpottarg,ifgradtarg,ifhesstarg);
%
% uu = U.pottarg.';
% uu = h*uu;
%
% end