function sol = solve_HH_eqn(sol,zzk, plot_panel)
% Solve_HH_eqn
% This code will solve the Helmholtz equation on a domain of choice. The
% outputs are intended to be used to solve the wave equation
%
% input is an object that specifies various experiment choices.
% THESE ARE DESCRIBED IN DETAIL IN THE TEST_DRIVER FILE.
%
% output is an object with various properties of interest:
%%%%%% INHERIT FROM INPUT OBJ:
% sol.bnd_type = boundary (scatterer type, string descrip, parameters)
% sol.inc = inc wave type (string descrip)
% sol.field_type = field type (total, scatter or inc)
% sol.solver_type = solver used for HH eqns
% sol.eval_type_freq = eval method used in freq space (e.g., 'fmm')
% sol.dom = boundary type (string descrip)
% sol.dom_range = range of domain.
% sol.dompts = spatial points for eval in domain(if supplied)
% sol.szIFT = # frequency pts needed in case of automatic sampling
% sol.NystM = disc. size for Nystrom systems.
% sol.domaingrid = {xx,yy,zz, keep}; (grid info if domain autosampled)
% sol.bandlimit = freql (lower and upper limit)
% sol.gam = boundary as a function handle
% sol.incF = incident field in freq. space

%%%%%%%%%%%%% COMPUTED HERE
% sol.freq = frequency points where solution has a value.
% sol.tau = density layer functions at each point in sol.freq
% sol.freqvals = Matrix of sol in freq dom at freqs
% sol.timsolv = timing for solve HH eqn
% sol.timevalfreq = timing for eval in freq space (not counting solve)
% sol.s = quadrature object from zeta-trap
%
% For now, the code is limited to incident waves for which we have
% explicit representations (function handles) in time domain.
%

%% Initialization
Uinc_pre = sol.incF;
M = sol.NystM;
s = sol.s;
W = sol.bandlimit;
solver_type = sol.solver_type;
eval_type_freq = sol.eval_type_freq;
eval_type_time = sol.eval_type_time;

if ~isa(Uinc_pre,'cell')
    Uinc = {Uinc_pre};
else
    Uinc = Uinc_pre;
end

%% Compute frequencies of interest.
field_type = 'scatter';
switch eval_type_time
    case 'freq only'
        ww = sol.eval_freqs;
        wts = 1;
        field_type = sol.field_type;
    case 'fastsinc'
        % we need equally spaced points on the frequency domain:
        % an odd number of points.
        N = sol.szIFT;W = sol.bandlimit;
        if ~mod(N,2)
            N = N+1;
        end
        ww = linspace(W(1), W(end), N+1); ww = ww(1:end-1).';
        wts = []; 
    case 'GLquad'
        % eval pts will be GL quad nodes.
        N = sol.szIFT;W = sol.bandlimit;
        [ww, wts] = legpts(N, W);
        if isempty(wts) %deal with single freq case
            wts = 1;
        end
    case 'complexify'

        ncomp = 50; %hard code # side cont. pts.
        % want to relate to size of del later.
        delt = sol.delt; 
        [wdel, wts] = legpts(ncomp, [0, delt]);%quadpts for side contours
        if isempty(wts) %deal with single freq case
            wts = 1;
        end
        N = sol.szIFT; W = sol.bandlimit;
        if ~mod(N,2)
            N = N+1;
        end
        wl = linspace(W(1), W(end), N+1); wl = wl(1:end-1).';
        ww = {W(2)+wdel*1i, wl+delt*1i, W(1) + wdel*1i};
end

sol.freq = ww;
sol.wts = wts;

%% Compute solution, for computed frequencies, to the Helmholtz density layer
%  and solution at given points zzk

time2solveHH = 0; time2evalfreqs = 0;

for i1 = 1:length(Uinc)
    switch eval_type_time
        case 'complexify'
            tim = tic;
            sol_HH_dens1 = solve_HH_density(Uinc{i1}, s, W, ww{1}, M, solver_type, plot_panel);
            sol_HH_dens2 = solve_HH_density(Uinc{i1}, s, W, ww{2}, M, solver_type, plot_panel);
            sol_HH_dens3 = solve_HH_density(Uinc{i1}, s, W, ww{3}, M, solver_type, plot_panel);
            timout = toc(tim);
            time2solveHH = time2solveHH + timout;
            
            sol.s = {sol_HH_dens1.s, sol_HH_dens2.s, sol_HH_dens3.s};
            tau{i1} = {sol_HH_dens1.tau,sol_HH_dens2.tau,sol_HH_dens3.tau};

            tim = tic;
            %eval_type_freq = 'basic';
            Usf1 = eval_sol_freqs(sol_HH_dens1, zzk, field_type, eval_type_freq, ww{1});
            Usf2 = eval_sol_freqs(sol_HH_dens2, zzk, field_type, eval_type_freq, ww{2});
            Usf3 = eval_sol_freqs(sol_HH_dens3, zzk, field_type, eval_type_freq, ww{3});
            timout = toc(tim);
            time2evalfreqs = time2evalfreqs + timout;

            Usf{i1} = {Usf1, Usf2, Usf3};
        otherwise
            % Get Density Layer
            tim = tic;
            sol_HH_dens = solve_HH_density(Uinc{i1}, s, W, ww, M, solver_type, plot_panel);
            timout = toc(tim);
            time2solveHH = time2solveHH  + timout;

            tau{i1} = {sol_HH_dens.tau};

            tim = tic;
            Usf{i1} = eval_sol_freqs(sol_HH_dens, zzk, field_type, eval_type_freq, ww);
            timout = toc(tim);
            time2evalfreqs = time2evalfreqs + timout;
    end
end

sol.tau = tau;
sol.freqvals = Usf;
sol.timsolv = time2solveHH;
sol.timevalfreq = time2evalfreqs;


end