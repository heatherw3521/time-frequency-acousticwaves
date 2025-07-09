function varargout = HH_density(ka, B,s, solve_type)
% find the density function tau for a given freq
% ka, gam = boundary (a chebfun), B = function of
% freq and spatial pts, should be -inc field on boundary.

ord = 32;
ieta = 1i*real(ka);
x = s.x;
if isa(B, 'function_handle')
    R = B(x,ka); %set up RHS
else
    R = B;
end
if size(R,1) == 1
    R = transpose(R);
end

dimension = 2;
if isfield(s, 'Zu')
    dimension = 3;
end
% check if right hand side is all zeros
if max(abs(R))<1e-16 & dimension==2
    varargout{1} = zeros(size(R));
    return
end

N = length(x);



%%
switch solve_type
    case 'basic'
        %use backslash w/out compression:
        if dimension == 2
            if abs(ka) > 1e-16
                Ad = HelmDLP(ka,s,s,[],[],ord);
                As = HelmSLP(ka,s,s,[],[],ord);
                A = eye(N)/2 + Ad - ieta * As;
            else
                Ad = LapDLPmat(s,s);      % D operator
                A = eye(N)/2+Ad + s.w.';  % 0.5+D
            end
        elseif dimension == 3
            if abs(ka) > 1e-16
                ord = 5; % ord 32 takes too much time and memory
                lptypes = {'s','d','sn'}; ZZ = Helm3dPatchZetaSparse_multi(ka,ord,lptypes,s);
                Zs = ZZ{1}; Zd = ZZ{2};
                Ad = Helm3dDLPmat(ka,s,s);
                As = Helm3dSLPmat(ka,s,s);
                A = eye(N)/2 + Ad - ieta * As + (Zd - ieta*Zs);
            else
                ord = 5; % ord 32 takes too much time and memory
                lptypes = {'d'}; ZZ = Helm3dPatchZetaSparse_multi(ka,ord,lptypes,s);
                Zd = ZZ{1};
                zs = [0.3,-0.9,0; 0.5,0.85,0; -0.99,-0.1,0].';
                nul = 1./vecnorm(s.x-zs(:,1))'; % place a source inside torus to cancel nullspace
                A = Lap3dDLPmat(s,s); A = A + Zd + 0.5*eye(size(A)) + s.w.*nul; 
            end
        end
        varargout{1} = A\R;
end
end





