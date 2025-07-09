function sol = solve_HH_density(F_inc, s, freql, bandsamples, M, solve_type, plot_panel)
%solve wave eqn with sound soft b.c., 
%F_inc = incident field (a function of freq and loc)
%gam = boundary of scatter (a chebfun)
% dom = square domain we want to solve on (scatter inside dom)
% bandsamples = freq samples, freq locs where we need to solve HH.
% freql = frequency band to solve HH eqn on : [W1, W2]; 
% M = discretization param for solving Nystrom eqn. 
%
% solution sol is returned as an object with following fields:
% sol.freq = frequency points where solution has a value.
% sol.tau = density layer functions at each point in sol.freq
% sol.bandlimit = freql (lower and upper limit)
% sol.incF = incident field in freq. space
% sol.gam = boundary 
% sol.s = quadrature object from zeta-trap
%
%  NOTE THAT AS OF NOW, SOLVER IS SIMPLE: O(M^3) for every pt in bandsamples. 
%% PART 1: set up spatial domain 
sol.bandlimit = freql; 
sol.incF = F_inc;
sol.gam = s.Z; 
%% PART 2: Solve for densities:
if isa(F_inc, 'function_handle')
    B = @(z,w) -F_inc(z,w); 
else
    B = -F_inc; %otherwise assume B is boundary data in a vector. 
end
 
%%%%%%%%%%%%%%%%%%OLD VERSION WE SIMPLY DO DISCRETIZATION FOR GL QUAD%%%%
%[w, wts] = legpts(N, freql); 
%if isempty(wts) %deal with single freq case
%    wts = 1; 
%end
%sol.freq = w; 
%sol.qw = wts; 
%%%%%%%%%%%%%%%%%%%%%%% NEW VERSION DISC over freqs is supplied
w = bandsamples; %rename to match legacy code. 
sol.freq = w; %store for object
N = length(w); 

% now we get the density functions for comb. field eqns. This is needed
% at each frequency: 
dimension = 2;
if isfield(s, 'Zu')
    dimension = 3;
end

if dimension == 2
    if length(M) == 1 
        s = setupquad(s, M);
    else
        s = ns_setupquad(s, M, plot_panel);
        M = length(s.x);
    end
    tau = zeros(M,N);
else
    s = quadr_doubleptr_patch2(s, [M,M]);
    tau = zeros(M^2,N);
end 
sol.s = s; 

switch solve_type
    case 'basic'
        if isa(B, 'function_handle')
            for j = 1:length(w)
                ka = w(j); 
                tau(:,j) = HH_density(ka, B,s, 'basic');
            end
        else
            for j = 1:length(w)
                ka = w(j); 
                tau(:,j) = HH_density(ka, B(:,j),s, 'basic');
            end 
        end
        sol.tau = tau; 
end

end