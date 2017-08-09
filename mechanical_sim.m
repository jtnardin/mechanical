clear all; clc

x = linspace(0,1,41);
dt = 1e-3;
t = 0:dt:.92;
dx = x(2)-x(1);

tn = length(t);
xn = length(x);

%interior
xint = 2:xn-1;
%boundary
xbd1 = [1];
xbdn = [xn];

%params
mu = 10;
E = 10;
lambda = 5;
k = .5;
D = 1e-3;
s = 1;

q = [mu;E;lambda;k;D;s];

%initial conditions

nIC = .5*(1-tanh(25*(x-.2)))';
rhoIC = .6-.5*.8*tanh(25*(x-.2))';
uIC = zeros(xn,1);
%vector of all entries
UIC = [nIC ; rhoIC ; uIC];


%define matrices for operations
%first order derivative

D1 = 1/(2*dx)*sparse([xint xint xbd1 xbd1 xbdn xbdn],...
    [xint+1 xint-1 xbd1 xbd1+1 xbdn xbdn-1],...
    [ones(1,xn-2) -ones(1,xn-2) -2 2 2 -2],xn,xn);

%second order derivative
D2 = 1/dx^2*sparse([xint xint xint],...
    [xint+1 xint xint-1],...
    [ones(1,xn-2) -2*ones(1,xn-2) ones(1,xn-2)],xn,xn); %%-2 2 -2 2

D2bd = sparse([1 xn],[1 xn],[1 1],xn,xn);

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

 

[t,U] = ode45(@(t,y) mechanical_ode(t,y,q,xn,dx,dt,D1,D2,D2bd,A_pos,A_pos_0,...
    A_pos_1,A_neg,A_neg_0,A_neg_1,xint),t,UIC);


v = diff(U(:,2*xn+1:end))/dt;

figure

for i = 1:99
    subplot(2,2,1)
    hold off
    plot(x,U(10*i+1,1:xn))
    hold on
    plot(x,v(10*i+1,:));
    title(num2str(i))
    
    subplot(2,2,2)
    hold off
    plot(x,U(10*i+1,xn+1:2*xn))
    hold on
    plot(x,v(10*i+1,:));
    
    subplot(2,2,3)
    plot(x,U(10*i+1,2*xn+1:3*xn))
    pause(.25)
end

