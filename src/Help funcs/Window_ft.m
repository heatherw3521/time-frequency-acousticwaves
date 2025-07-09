function Res = Window_ft(f,zz,ww,s,N,H,varargin)
%Function you want to Transform: f
%Sectioning: s
%Frequency: w or Time: t (depending on output_type)

% Output F_k^slow where each row corresponds to windowd ft of around s_n and collumn to w
Res = zeros(length(zz),length(ww));
if size(zz,1) == 3 % We need something better for this
    nn = size(zz,2);
else
    nn = length(zz);
end
for i1 = 1:nn
    if size(zz,1) == 3 % We need something better for this
        z = zz(:,i1);
    else
        z = zz(i1);
    end
    for i2 = 1:length(ww)
        w = ww(i2);
        if length(s) == 1
            fk = f_windowed(@(t)f(z,t),H,s);
            if ~isempty(varargin)
                Res(i1) = Res(i1) + fk(w+s);
            else
                ff = @(t) fk(t+s) .* exp(1i*w.*t);

                I = [-H,H];

                Res(i1,i2) = Clen_Curt(ff,I,N);
            end
        end
    end
end
end

function fk = f_windowed(f,H,s)
fk = @(t) f(t).*window(t-s,H);
end

function w = window(tt,H)
w = tt * 0;
eta = @(u) exp(2*exp(-1./u)./(u-1));
for i1 = 1:length(tt)
    t = tt(i1);
    if (-H< t && t< -H/2)
        w(i1) = 1-eta(2*(t+H)/H);
    elseif (-H/2<= t && t<= H/2)
        w(i1) = 1 ;
    elseif (H/2<= t && t< H)
        w(i1) = eta((2*t-H)/H);
    end
end
end
