clear all; clc

x = linspace(-10,1,201);
dt = 1e-3;
t = 0:dt:4;
dx = x(2)-x(1);

tn = length(t);
xn = length(x);

%interior
xint = 2:xn-1;
%boundary
xbd1 = [1];
xbdn = [xn];
% % 
% % %params
% % mu = 1;
% % E = 1;
% % lambda = 1;
% % k = 1;
% % D = 1;
% % s = 1;
% 
% q = [mu;E;lambda;k;D;s];


%%%% Using this to model base model from Tranquillo-Murray 1992. Something
%%%% still wrong with my boundary conditions, apparently. When s small,
%%%% results look somewhat similar (but s = 1 in their paper..)
r0 = 1;
mu = 1;
tau = 1/2;
gamma = 1;
s = .001;

q = [r0, mu, tau, gamma, s];

%initial conditions

nIC = .5*(1-tanh(25*(x-.2)))';
rhoIC = 1*ones(length(x),1);%.6-.5*.8*tanh(25*(x-.2))';
uIC = zeros(xn,1);
%vector of all entries
UIC = [nIC ; rhoIC ; uIC];


%define matrices for operations
%first order derivative

D1 = 1/(2*dx)*sparse([xint xint xbd1 xbd1 xbdn xbdn],...
    [xint+1 xint-1 xbd1 xbd1+1 xbdn xbdn-1],...
    [ones(1,xn-2) -ones(1,xn-2) 0 0 0 0],xn,xn);

%second order derivative
D2 = 1/dx^2*sparse([xint xint xint],...
    [xint+1 xint xint-1],...
    [ones(1,xn-2) -2*ones(1,xn-2) ones(1,xn-2)],xn,xn); %%-2 2 -2 2

D2bd_D = sparse([1 xn],[1 xn],[1 1],xn,xn);

D2bd_N = sparse([1 1 xn xn],[1 2 xn-1 xn],[-2 2 2 -2],xn,xn);

% Matrices for advection with flux limiters

A_pos = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[vw.*(-1+sw/2); ...
    (ve.*(1-1*se/2)-vw.*sw/2); (ve.*se/2)],xn,xn);

A_pos_0 = @(se,ve,ind,dn) sparse([ind ind],[ind ind+dn],[ve.*(1-1*se/2); ve.*se/2],xn,xn);

A_pos_1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[vw.*(-1+sw/2); ...
    (-vw.*sw/2)],xn,xn);

A_neg = @(se,sw,ve,vw,ind,dn) sparse([ind ind ind],[ind-dn ind ind+dn],[(-vw.*sw/2); ...
    (ve.*se/2+vw.*sw/2-vw); (ve-ve.*se/2)],xn,xn);

A_neg_1 = @(sw,vw,ind,dn) sparse([ind ind],[ind-dn ind],[(-vw.*sw/2); ...
    (vw.*sw/2-vw)],xn,xn);

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