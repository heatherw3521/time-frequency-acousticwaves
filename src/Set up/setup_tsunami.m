function [uinc, Uinc] = setup_tsunami(sol,T, w0, sigma,c,r,W_est, dim, zzk)
% parameters to set up the Gaussian packet incident wave
% 
% T = shift in time
% w0 = shift in phase
% sigma = width parameter (smaller sigma = bigger band in time/smaller band in freq)
% c = wavespeed
% r = dir. of propogation (number in complex plane assoc. with vector in
% R2)
% W_est = [a, b]: estimated band on which |Uinc| is significant
%%
uinc_part = @(t) 1/(sqrt(2*pi)*sigma) * exp(-1i*w0*t) .* exp(-t.^2/(2*sigma^2));
Uinc_part = @(w,A) exp(-(w-w0).^2*sigma^2/2) .* exp(1i*A.*w);
if dim==2
    r = r/abs(r); %incident field direction of prop.
    Uinc =  @(z, w) Uinc_part(w,real(r).*real(z)/c + imag(r).*imag(z)/c+T);
    uinc = @(z, t) uinc_part(t-(real(r).*real(z)/c + imag(r).*imag(z)/c+T));
%% check if the function is small off the freq band: assumes w_est(1) > 0
%
a = W_est(1); b = W_est(2); 
happy = false; num1 = 0; 
rp = -2+1i; % hopefully outside the scatterer. Eventually replace 
% with a way to fix a test point in domain.
while ~happy
    num1 = num1+1; 
    tst = linspace(a/2, a, 200); 
    if max(abs(Uinc(tst,rp))) < 1e-11 || num1 >= 20
        happy = true; 
    else
        a = a/2; 
    end
end
W(1) = a; 
happy = false; 
num2 = 0; 
while ~happy
    num2 = num2+1; 
    tst = linspace(b, (b+ 5)/2, 200); 
    if max(abs(Uinc(tst,rp))) < 1e-11 || num2 >= 20
        happy = true; 
    else
        b = (b+5)/2; 
    end
end
W(2) = b; 
%
if max(num1, num2) > 19
    warning('Incident wave support is larger than the frequency band')
end
%%
% check if IC is satisfied: Is uinc zero valued on all domain at t =  0?
% also check if uinc passes through domain and everything is zero-valued at
% 2*T (the final value of the timeband). 
zpts = sol.domaingrid{3};
ic = max(abs(uinc(zpts,0))); 
if ic > 1e-10
    warning('initial conditions not satisfied')
end

num1 = 0; 
happy = false; 
while ~happy
    fc = max(abs(uinc(zpts, 2*T))); 
    if fc < 1e-10 || num1 == 5
        happy = true; 
    else
        T = T + 10;
    end
end
if num1 > 4
    warning('incident wave does not complete domain passthrough in given time')
end

% TO DO: 3D version (for now just sets up the basic thing)
elseif dim==3
     r = r/sqrt(sum(r.^2)); %incident field direction of prop.
     Uinc = @(z, w) Uinc_part(w,r(1).*z(1,:)/c + r(2).*z(2,:)/c + r(3).*z(3,:)/c+T);
     uinc = @(z, t) uinc_part(t-(r(1).*z(1,:)/c + r(2).*z(2,:)/c + r(3).*z(3,:)/c+T));
end


end