function u = eval_sol_time(w,W,zz, Uf, t, evalstyle, delta, field_type, quad, inc_time)
%Evaluates the Fourier transform of the BB HH solution at the time T: 
%
% if sol(zz,w) denotes the HH solution at freq w, location zz, then the 
% computed solution to the wave equation in time is given by its inverse
% Fourier transform: 
%
% u(zz,t) = 1/2/pi * ( int_{bandlimit} sol(w, zz) * exp(-1i*w*t) dw. 
%
% Uf = values of sol for all pts (w,zz). Call eval_sol_freqs_faster
% to generate this. 
%
% field_type = 'scatter', 'total', or 'incident' 
%
% we provide several ways to do the evaluation; the way you select should 
% correspond to how the frequency points in sol were sampled:
%
% evalstyle = 'fastsinc'. This option constructs a trigonometric polynomial
% approximation to sol(w,zz) over w using the fft. This leads to a closed form
% solution to the inverse involving sums of sinc functions. This is the
% default choice and should be used if the scatterer does not have
% cavities/trapping qualities.
%
% evalstyle = 'complexify'. This option perturbs w into the complex plane by a
% user-supplied factor delta, and sets up a rectangular contour inside
% of which sol is analytic. We use the fast sinc transform to evaluate the 
% top integral, and orthogonalization to handle the side contours, which 
% grow exponentially. 
%
% evalstyle = 'AAA'. This option uses AAA to construct a trigonometric
% rational approximation to the sol object. The points used to 
% run AAA will coincide with those in sol.freqs. The rational approximation 
% allows us to write a solution to the integral in closed form as a
% complex-valued exponential sum that must be evaluated.
%
%
% evalstyle = 'GLquadrature'. This option simply applies a quadrature rule to
% the integral to evaluate. The user should supply quadrature weights. The
% quadrature points are assumed to coincide with those in sol.freq. 

%%
% quad = quadrature weights. quad = [] if you are not using quadrature
%
% inc_time = function handle u_inc(zz,t) of incident wave in time domain. 
% only needed if field_type = 'total'. Otherwise, can set to []. 
%%
if ~(size(zz,1) == 3)
    zz = zz(:);
    nn = length(zz);
else
    nn = size(zz,2); 
end

t = t(:); t = t.';

if strcmpi(field_type, 'incident')
    u = zeros(length(zz),length(t)); % Dumb solution, needs to improve in definition of incident wave
    for i1 = 1:length(t)
        u(:,i1) = inc_time(zz,t(i1));
    end
    return
end

switch evalstyle
    case 'fastsinc'
        % use fast sinc transform to evaluate the inv. Fourier transform. 
        w = w(:); 
        M = length(w); 
        W1 = W(1);
        P = W(2) - W(1); 
        % if M is even, we need to fix it: 
        % for now pad with zero: 
        if ~mod(M,2)
            w = [w;0];
            M = length(w);
        end
        %u = zeros(length(zz),length(t)); %solution store
        %%
        % apply FFT to get trig coeffs for Uf(w,\cdot): 
        %step 1: get fft coeffs for f
        mm = (M-1)/2;  
        C =fft(Uf, [], 2)/M; C = fftshift(C,2); 
        % we now want to compute the following at each zz:
        % G_k = e^(-1i*t/2)*sum_{j = -mm to mm} q_j sinc(t_k/2/pi-j)
        TT = fastsincwrapper(mm, t, P, C, nn);
        u = repmat(exp(-1i*(P/2 + W1)*t),nn,1).*TT.';
        u = u*(P/(2*pi)); 

    case 'GLquad'
        wts = quad; 
        wts = wts(:); 
        u = zeros(nn, length(t)); %rows vary by location, cols vary by time. 
%do inv_fourier transform to get time information on scattered field
        for j = 1:length(t)
            tt = t(j); 
            u(:,j)= (Uf.*exp(-1i*w.'*tt))*wts;
        end
        u = u/(2*pi);
        % to do next
 %%       
    case 'complexify'
     % we evaluate the solution on a complex contour: 
     % [w2, w1] + delta*1i 
     %
     % get the freq pts
     wc2 = w;  %ordered from w1 to w2
     W1 = W(1); W2 = W(2);
     P = W2-W1;
      
     M = length(wc2);
     % if M is even, we need to fix it: 
     % for now pad with zero: 
     if ~mod(M,2)
        wc2 = [wc2;0];
        M = length(wc2);
     end
     mm = (M-1)/2;  
     C = fft(Uf, [], 2)/M; C = fftshift(C,2); 
     % we now want to compute the following at each zz:
     % G_k = e^(-1i*t/2)*sum_{j = -mm to mm} q_j sinc(t_k/2/pi-j)
     TT = fastsincwrapper(mm, t, P, C, nn);
     u = repmat(exp(-1i*(P/2 + W1 + 1i*delta)*t),nn,1).*TT.';
     u = -u*(P/(2*pi)); 

     % %now we do c1 and c3 via GLquad:
     % %quad = N;%length of quadrature rule
     % %[pts, wts] = legpts(N,0, delta);
     % wts = quad; wts = wts(:); 
     % for j = 1:length(t)
     %     tt = t(j); 
     %     u1(:,j)= (Ufc1.*exp(wc1.'*tt))*wts;
     %     u3(:,j)= (Ufc3.*exp(wc3.'*tt))*wts; 
     % end
     % u1 = u1/(2*pi);
     % u3 = -u3/(2*pi);
     % % now form u1: 
     % u = -(u1 + u2 + u3);
end

% if strcmpi(field_type, 'total')
%     if size(zz,1) == 3
%         for i1 = 1:length(t)
%             u(:,i1) = u(:,i1) + transpose(inc_time(zz,t(i1))); 
%         end
%     else
%         u = u + inc_time(zz,t); 
%     end
% end
end
%%%%%%%%%%%%%%%%%%%%%%%subfunctions%%%%%%%%%%%%%%%%%
function TT = fastsincwrapper(mm, t, P, C, nn)
        ifl = 1; 
        a1 = t.'*(P/(2*pi)); a1 = a1(:); 
        klocs = (-mm:mm);%klocs = klocs(:); 
        Q = C.*repmat((-1).^(klocs),nn,1); Q = Q.';  
        tol = 1e-12; 
        mode = 'legendre';
        %for j = 1:nn
        %    u(j,:) = exp(-1i*t.'/2).*sinc1d(ifl,a1,klocs,Q(:,j),tol,mode);
        %end
        %u = u/(2*pi); 
        TT = sinc1d(ifl,a1,klocs,Q,tol,mode);
end
