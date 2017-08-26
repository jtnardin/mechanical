clear all; clc

%define grids
xn = 21;
yn = 21;
total = xn*yn;

x = linspace(-1,1,xn);
y = linspace(0,1,yn);

dx = x(2)-x(1);
dy = y(2)-y(1);
dt = 1e-3;
t = 0:dt:1;
tn = length(t);
[Y,X] = meshgrid(y,x);


%construct boundary points, interior
y_bd_0 = 1:xn;
y_bd_l = xn*(yn-1)+1:total;

y_bd = union(y_bd_0,y_bd_l);

x_bd_0 = 1:xn:total;
x_bd_l = xn:xn:total;

x_bd = union(x_bd_0,x_bd_l);

bd = union(y_bd,x_bd);

%interior for x
x_int = 1:total;
x_int(x_bd) = [];

%interior for y
y_int = 1:total;
y_int(y_bd) = [];

%interior for all
xy_int = 1:yn*xn;
xy_int(bd) = [];


%%%% Using this to model base model from Tranquillo-Murray 1992. Something
%%%% still wrong with my boundary conditions, apparently. When s small,
%%%% results look somewhat similar (but s = 1 in their paper..)
r0 = 1;
mu1 = 1;
mu2 = 1;
tau = 1/2;
E = 1;
nu = .25;
gamma = 1;
s = 1;
D = 0.1;

q = [mu1,mu2,E,nu,s,tau,gamma,r0,D];

%initial conditions

nIC = .5*(1-tanh(25*(X-.2)));
nIC = nIC(:);
rhoIC = 1*ones(size(X));%.6-.5*.8*tanh(25*(x-.2))';
rhoIC = rhoIC(:);
uIC = zeros(size(X));
uIC = uIC(:);
%vector of all entries
UIC = [nIC ; rhoIC ; repmat(uIC,2,1)];


%define matrices for operations
%first order derivative

D1X = @(ind,n) 1/(2*dx)*sparse([ind ind],...
    [ind+1 ind-1],[ones(1,length(ind)) -ones(1,length(ind))],n,n);

D1Y = @(ind,n) 1/(2*dy)*sparse([ind ind],...
    [ind+xn ind-xn],[ones(1,length(ind)) -ones(1,length(ind))],n,n);

%%%%need something for boundary (currently using 0 at all boundary points
%%%%in 1d)

%second order derivative for velocity computation
D2X = @(ind,n) 1/dx^2*sparse([ind ind ind],...
    [ind+1 ind ind-1],...
    [ones(1,length(ind)) -2*ones(1,length(ind)) ones(1,length(ind))],n,n); %%-2 2 -2 2

D2Y = @(ind,n) 1/dy^2*sparse([ind ind ind],[ind+xn ind ind-xn],...
    [ones(1,length(ind)) -2*ones(1,length(ind)) ones(1,length(ind))],n,n); %%-2 2 -2 2

D2XY = @(ind1,ind2,n) 1/(4*dx*dy)*sparse([ind1 ind1 ind1 ind1],...
    [ind2+xn+1 ind2+1-xn ind2-1+xn ind2-1-xn],...
    [ones(1,length(ind1)) -1*ones(1,2*length(ind1)) ones(1,length(ind1))],n,n); %%-2 2 -2 2


D2bd_D = @(ind,n) sparse(ind,ind,ones(length(ind),1),n,n);


D2Xbd_N = @(ind0,ind1,n) sparse([ind0 ind0 ind1 ind1],[ind0 ind0+1 ind1-1 ind1]...
    ,[-2*ones(1,length(ind0)) 2*ones(1,length(ind0)) 2*ones(1,length(ind1))...
    -2*ones(1,length(ind1))],n,n);


D2Ybd_N = @(ind0,ind1,n) sparse([ind0 ind0 ind1 ind1],[ind0 ind0+xn ind1-xn ind1]...
    ,[-2*ones(1,length(ind0)) 2*ones(1,length(ind0)) 2*ones(1,length(ind1))...
    -2*ones(1,length(ind1))],n,n);

% Matrices for advection with flux limiters

A_pos = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[vw.*(-1+sw/2); ...
    (ve.*(1-1*se/2)-vw.*sw/2); (ve.*se/2)],total,total);

A_pos_0 = @(se,ve,ind,dn) sparse([ind ind],[ind ind+dn],[ve.*(1-1*se/2); ve.*se/2],total,total);

A_pos_1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[vw.*(-1+sw/2); ...
    (-vw.*sw/2)],total,total);

A_neg = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-vw.*sw/2); ...
    (ve.*se/2+vw.*sw/2-vw); (ve-ve.*se/2)],total,total);

A_neg_1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[(-vw.*sw/2); ...
    (vw.*sw/2-vw)],total,total);

A_neg_0 = @(se,ve,ind,dn) sparse([ind ind],[ind ind+dn],[ve.*se/2; (ve-ve.*se/2)],xn,xn);



tic

[t,U] = ode15s(@(t,y) mechanical_ode_TM_base_2d(t,y,q,dx,dy,D1X,D1Y,D2X,D2Y,D2XY,...
    D2bd_D,D2Xbd_N,D2Ybd_N,A_pos,A_pos_0,A_pos_1,A_neg,A_neg_0,A_neg_1,xy_int,bd,x_bd_0,x_bd_l,...
    y_bd_0,y_bd_l,total,x_int,y_int,xn,yn),t,UIC);

toc

tn = length(t);


%%%%%%%% plot

%extract each variable to its matrix (t,x,y)
n = reshape(U(:,1:total),tn,xn,yn);
rho = reshape(U(:,total+1:2*total),tn,xn,yn);
ux = reshape(U(:,2*total+1:3*total),tn,xn,yn);
uy = reshape(U(:,3*total+1:4*total),tn,xn,yn);
vx = diff(ux,1)/dt;
vy = diff(uy,1)/dt;


figure('units','normalized','outerposition',[0 0 1 1])

for i = 1:10:tn-1
   subplot(2,2,1)

   
   hold off
   contourf(Y,X,squeeze(n(i,:,:)),'edgecolor','none')
%    view(2)
   xlabel('y')
   ylabel('x')
   title(['cells, t = ' num2str(t(i))])
   
   caxis([0 1])
   colorbar

   hold on
   
   quiver(Y,X,squeeze(vy(i,:,:)),squeeze(vx(i,:,:)),'linewidth',2,'color','r');
   
   subplot(2,2,2)
   
   contourf(Y,X,squeeze(rho(i,:,:)),'edgecolor','none')
   view(2)
   xlabel('y')
   ylabel('x')
   title(['collagen, t = ' num2str(t(i))])
   
   
%    caxis([min(min(min(rho))) max(max(max(rho)))]);
   colorbar
   
   subplot(2,2,3)
   
   contourf(Y,X,squeeze(ux(i,:,:)),'edgecolor','none')
   view(2)
   xlabel('y')
   ylabel('x')
   title(['ux, t = ' num2str(t(i))])
   
   
   caxis([min(min(min(ux))) max(max(max(ux)))]);
   colorbar
   
   subplot(2,2,4)
   
   contourf(Y,X,squeeze(uy(i,:,:)),'edgecolor','none')
   view(2)
   xlabel('y')
   ylabel('x')
   title(['uy, t = ' num2str(t(i))])
   
   
   caxis([min(min(min(uy))) max(max(max(uy)))]);
   colorbar
   
   pause(.125)
    
end
