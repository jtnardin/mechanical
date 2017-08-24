clear all; clc

%define grids
xn = 41;
yn = 11;
total = xn*yn;

x = linspace(-1,1,xn);
y = linspace(0,1,yn);

dx = x(2)-x(1);
dy = y(2)-y(1);
dt = 1e-3;
t = 0:dt:4;
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


%extra bound for sensors -- sensing second, second to last points from
%INTERIOR points.
y_bd_1_int = y_bd_0 + xn;
y_bd_nm1_int = y_bd_l - xn;

x_bd_1_int = y_bd_0 + 1;
x_bd_nm1_int = x_bd_l - 1;




%%%% Using this to model base model from Tranquillo-Murray 1992. Something
%%%% still wrong with my boundary conditions, apparently. When s small,
%%%% results look somewhat similar (but s = 1 in their paper..)
r0 = 1;
mu1 = 1;
mu2 = 1;
tau = 1/2;
E = 1;
nu = .75;
gamma = 1;
s = .01;

q = [mu1,mu2,E,nu,s,tau,gamma];

%initial conditions

nIC = .5*(1-tanh(25*(X-.2)));
nIC = nIC(:);
rhoIC = 1*ones(size(X));%.6-.5*.8*tanh(25*(x-.2))';
rhoIC = rhoIC(:);
uIC = zeros(size(X));
uIC = uIC(:);
%vector of all entries
UIC = [nIC ; rhoIC ; uIC];


%define matrices for operations
%first order derivative

D1X = @(ind) 1/(2*dx)*sparse([ind ind],...
    [ind+1 ind-1],[ones(1,length(ind)) -ones(1,length(ind))],2*total,2*total);

D1Y = @(ind) 1/(2*dy)*sparse([ind ind],...
    [ind+xn ind-xn],[ones(1,length(ind)) -ones(1,length(ind))],2*total,2*total);

%%%%need something for boundary (currently using 0 at all boundary points
%%%%in 1d)

%second order derivative for velocity computation
D2X = @(ind) 1/dx^2*sparse([ind ind ind],...
    [ind+1 ind ind-1],...
    [ones(1,length(ind)) -2*ones(1,length(ind)) ones(1,length(ind))],2*total,2*total); %%-2 2 -2 2

D2Y = @(ind) 1/dy^2*sparse([ind ind ind],[ind+xn ind ind-xn],...
    [ones(1,length(ind)) -2*ones(1,length(ind)) ones(1,length(ind))],2*total,2*total); %%-2 2 -2 2

D2XY = @(ind1,ind2) 1/(4*dx*dy)*sparse([ind1 ind1 ind1 ind1],...
    [ind2+xn+1 ind2+1-xn ind2-1+xn ind2-1-xn],...
    [ones(1,length(ind1)) -1*ones(1,2*length(ind1)) ones(1,length(ind1))],2*total,2*total); %%-2 2 -2 2


D2bd_D = @(ind) sparse(ind,ind,ones(length(ind),1),2*total,2*total);


v = velocity_solve(q,D1X,D1Y,D2X,D2Y,D2XY,D2bd_D,rhoIC,nIC,repmat(uIC,2,1)...
    ,xy_int,bd,total);


D2bd_N = sparse([1 1 xn xn],[1 2 xn-1 xn],[-2 2 2 -2],total,total);

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

[t,U] = ode15s(@(t,y) mechanical_ode_TM_base_w_BC(t,y,q,xn,dx,dt,D1,D2,D2bd_D,D2bd_N,A_pos,A_pos_0,...
    A_pos_1,A_neg,A_neg_0,A_neg_1,xint),t,UIC);

toc

v = diff(U(:,2*xn+1:end))/dt;

%video
% %     
% figure
% hold on
% 
% dt_plot = floor(length(t)/200);
% 
% for i = 1:199
% %     subplot(2,2,1)
%     hold off
%     plot(x,U(i*dt_plot,1:xn))
%     hold on
%     plot(x,v(i*dt_plot,:));
%     title(num2str(t(i*dt_plot)))
%     
% %     subplot(2,2,2)
% %     hold off
%     plot(x,U(i*dt_plot,xn+1:2*xn))
% %     hold on
% %     plot(x,v(10*i+1,:));
%     
% %     subplot(2,2,3)
%     plot(x,U(i*dt_plot,2*xn+1:3*xn))
%     
%     axis([x(1) 1 -.3 1.3])
%     pause(.25)
%     
%     
% end


%subplot
figure('unit','normalized','outerposition',[0 0 1 1])

dt_plot = floor(length(t)/4);

subplot(2,2,1)
hold on

for i = 1:dt_plot:length(t)
    plot(x,U(i,1:xn),'b')
end

xlabel('x')

axis([x(1) 1 -.3 1.3])

subplot(2,2,2)
hold on

for i = 1:dt_plot:length(t)
    plot(x,U(i,xn+1:2*xn),'b')
end

xlabel('x')

axis([x(1) 1 .9 1.1])

subplot(2,2,3)
hold on

for i = 1:dt_plot:length(t)
    plot(x,U(i,2*xn+1:3*xn),'b')
end

xlabel('x')

axis([x(1) 1 -2 0])