function f = freq2time_fastsinc(F,xx,t,Ib, m, tol)
% This code sets up a function
% % handle that uses the fast sinc transform to evaluate
% a potentially oscillatory integral of the form
%
% f(x, t) = 1/2/pi * \int_{a}^{b} F(x,w) e^{-iwt} dw.
%
% Here, x is a point in space, t = time, w= frequencies. This is the 
% transform from the frequency domain to the time domain. 
%
% For now, we assume that F(x,w) is periodic (choose a, b so that
% F(x,w) is approx 0 off the interval [a, b]). 
%
% input: F = function handle F(x,w), x = spatial pts, t = time pts
% Ib = integration bounds, m = approximate bandlimit of F(x,w) 
% (see memo_wv_2), tol = tolerance parameter (rough idea of rel. acc)
%  
%
% output: function values f(j) = f(x, t(j)).
%
% TODO: make work for multiple x vals.  


%%
% compute Fourier coeffs for F(x,w) rescaled to live on [0, 1). 
% we assume for a fixed x, F(x,w) \approx \sum_{j = -m}^m c_j exp(2*pi*1i*j*w);
%TO DO, set up for many x locs. 

a = Ib(1); b = Ib(2); 
P = (b-a);
ww = @(y) P*y + a; 

%loop tries to make sure we resolve F(x,w) with enough Fourier coeffs.
happy = 0; ct = 0; ctmax = 4; 
while ~happy
    ct = ct +1;
    M = 2*m + 1; %total number of Fourier coeffs in the series 
    pts = a + P*(0:M-1)'/M; 
    c =fft(F(xx,pts))/M; 
    %cutoff = StandardChop(c(1:m),tol);
    cutoff = standardChop(c(1:m),tol);
    if cutoff < m 
        happy =1;
    elseif ct > ctmax
        warning("we are not happy.")
        break
    else
        m = 2*m; 
    end
end
c = fftshift(c); 

%%
% the integral f(t,x) can be expressed as follows (see memo wv_2): 
%
% P*exp(it*(a-P/2)) * sum_{j = -m}^{m} c(j)*(-1)^j * sinc(P*t/2/pi - j), 
%
% where sinc(x) = sin(pi*x)/pi/x. 
%
% To evaluate the sum, we use the fast sinc transform: 
ifl = 1; 
a1 = P*t/2/pi; a1 = a1(:); 
klocs = (-m:m); 
q = c.*(-1).^(klocs'); 
tol = 1e-12; 
mode = 'legendre';
f = P*exp(1i*t*(a-P/2)).*sinc1d(ifl,a1,klocs,q,tol,mode);
f = f/2/pi; 
end











