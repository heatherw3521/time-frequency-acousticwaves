function res = Clen_Curt(gg,I,n,varargin)

if isempty(varargin)
    g = gg;
    res = 0;
else
    zz = varargin{1};
    if (size(zz,1) == 1 || size(zz,2) == 1) && ~isreal(zz(1))
        z = zz(2:end);
            g = @(t) gg(zz(1),t);
    res = zeros(length(zz),1);
    else
        z = zz(:,2:end);
        g = @(t) gg(zz(:,1),t);
        res = zeros(size(zz,2),1);
    end
end
if n == 0
    a = I(1);
    b = I(2);
    
    mapp = @(x) (b+a)/2 + x * (b-a)/2;
    Dmapp =  (b-a)/2;
    f = @(x) g(mapp(x)) * Dmapp;
    
    fc = UltraFun(0,f,500, 1e-13);
    
    n = length(fc);
    w = zeros(n,1);
    for k = 0:(n-1)
        if mod(k,2) == 0
            w(k+1) = 2. / (1. - k^2) * sqrt(2.);
        end
    end
    w(1) = 2;
    w(end-1:end) = w(end-1:end)/2;
    res(1) = w' * fc;
    if ~isempty(varargin)
        for i2 = 2:size(zz,2)
            g = @(t) gg(zz(:,i2),t);
            f = @(x) g(mapp(x)) * Dmapp;
            fc = UltraFun(0,f,500, 1e-13);
            
            n = length(fc);
            w = zeros(n,1);
            for k = 0:(n-1)
                if mod(k,2) == 0
                    w(k+1) = 2. / (1. - k^2) * sqrt(2.);
                end
            end
            w(1) = 2;
            w(end-1:end) = w(end-1:end)/2;
            res(i2) = w' * fc;
        end
    end
else
    if isempty(varargin)
        res = 0;
        g = gg;
    else
        z = varargin{1};
        res = zeros(length(z),1);
        g = @(t) gg(z,t);
    end

    a = I(1);
    b = I(2);
    N=n-1; bma=b-a;
    c=zeros(n,2);
    c(1:2:n,1)=(2./[1 1-(2:2:N).^2 ])'; c(2,2)=1;
    f=real(ifft([c(1:n,:);c(N:-1:2,:)]));
    w=bma*([f(1,1); 2*f(2:N,1); f(n,1)])/2;
    x=0.5*((b+a)+N*bma*f(1:n,2));

    for i1 = 1:length(x)
        func_eval = g(x(i1));
        if size(func_eval,2) ~= 1
            func_eval = transpose(func_eval);
        end
        res = res + func_eval .* w(i1);
    end
end
end

function [avec, bvec] = Jacobi_ab(a,b,N)
NS = 0:(N-1);
avec = zeros(N,1); bvec = zeros(N,1);
for i1 = NS+1
    n = NS(i1);
    if a+b==-1 && n==0
        bvec(i1) = sqrt(2*a*b);
    else
        bvec(i1) =2*sqrt(n+1).*sqrt(n+a+1).*sqrt(n+b+1).*sqrt(n+a+b+1)./((2*n + a +b + 2).*sqrt(2*n + a +b + 3).*sqrt(2*n + a + b + 1));
    end
    if (a+b==0 || a+b==-1) && n==0
        avec(i1) =  (b-a)/(a+b+2);
    else
        avec(i1) = (b^2 - a^2)./((2*n + a +b + 2).*(2*n + a +b));
    end
end
end

function fcoeff = UltraFun(lambda,f,N, tol)
n = 8;
while N>n
    [a,b] = Jacobi_ab(lambda-.5,lambda-.5,n);
    b = [0;b];
    E = spdiags([b(2:end),a,b(1:n)],-1:1,length(a),length(a));
    [Evec, Eval] = eigs(E,n);
    Eval = diag(Eval);
    fcoeff = Evec*(diag(Evec(1,:))*f(Eval));
    if sum(abs(fcoeff(end-4:end))) < tol
        return
    end
    n = n*2;
end
warning(['Maximal DOF reached, precision: ', num2str(sum(abs(fcoeff(end-4:end))))])
end

