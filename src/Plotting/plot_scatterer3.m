function plot_scatterer3(gam, ht, fillit)
% plot a scatterer on a 3D plot at the height ht
% input should be a chebfun gam defined by a parametrization t \in [0, 2pi]

% here we fatten a bit to cover the raggedy cutoff of the mesh.  
% this is bespoke and not a final fix to the problem
if isa(gam, 'chebfun')
    n = normal(gam);
    %n = @(x) ([x,x]).^0-1;
else
    n = @(x) ([x,x]).^0-1;
end
tpts = linspace(0,2*pi, 8001).'; tpts = tpts(1:end-1);
bpts = gam(tpts); nn = n(tpts); nn = nn./sqrt(nn(:,1).^2 + nn(:,2).^2); nn(isnan(nn)) = 0;
%nn = nn(:,1) + 1i*nn(:,2); gbts = bpts + .08*nn; 

nn = nn(:,1) + 1i*nn(:,2); gbts = bpts + .02*nn;

bptsx = real(gbts);
bptsy = imag(gbts);
%bptsx = real(gam(bpts));
%bptsy = imag(gam(bpts));
mydisc_plot(bptsx.',bptsy.', ht)
%%
% add gray: 
%
if fillit ==1
    h = mydisc_fill(bptsx, bptsy,[1 1 1], ht);
    %h.FaceAlpha = .75; %some transparency
end



end


function mydisc_plot(x,y, ht)
n = length(x);
dz = sqrt((x(2:n) - x(1:(n-1))).^2+(y(2:n) - y(1:(n-1))).^2);
ex_I = find(dz>=.5);
plt_points = [0;ex_I(:);n];
hold on
for i1 = 1:(length(plt_points)-1)
    indi = [(plt_points(i1)+1):plt_points(i1+1),(plt_points(i1)+1)];
    plot3(x(indi), y(indi), ht*ones(length(indi),1), '-k', 'linewidth', 4);
end
%Problem for circles
%xlim([-5, 5])% I don't think we want to set this 
%ylim([-5,5])
%
end

function h = mydisc_fill(x, y, clr, ht)
n = length(x);
dz = sqrt((x(2:n) - x(1:(n-1))).^2+(y(2:n) - y(1:(n-1))).^2);

ex_I = find(dz>=.5);
plt_points = [0;ex_I(:);n];

hold on
for i1 = 1:(length(plt_points)-1)
    indi = [(plt_points(i1)+1):plt_points(i1+1),(plt_points(i1)+1)];
    h = fill3(x(indi), y(indi),ht*(ones(size(indi))),clr);
    %h.FaceAlpha = .75; %some transparency
end
end