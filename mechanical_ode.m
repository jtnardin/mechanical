function yprime = mechanical_ode(t,y,q,xn,dx,dt,D1,D2,D2bd,A_pos,A_pos_0,...
    A_pos_1,A_neg,A_neg_0,A_neg_1,xint)

%     t

    n = y(1:xn);
    rho = y(xn+1:2*xn);
    u = y(2*xn+1:3*xn);
    
    
    %solve for v, which also gives us du/dt
    v = -1/q(1)*(D2 + D2bd)\(q(2)*((D1*rho).*(D1*u) + (rho).*(D2*u)) + ...
        q(3)*((D1*rho).*n + rho.*(D1*n)) - q(6)*rho.*n);
%     

    %where is velocity positive, negative
    
    v_pos = v(2:end-1)>=0;
    v_neg = v(2:end-1)<0;
    
    x_int_p = xint(v_pos);
    x_int_n = xint(v_neg);
    
    %velocities for computation
    ve = v(2:end);
    ve0 = v(1);
    vw = v(1:end-1);
    vw_n = v(end);

    %implicitness of scheme
    theta = 0.5;


    %sigma for flux limiters
    sigma = @(r) (r+abs(r))./(1+abs(r));

    
    %construct Apos (positive matrix)
    
    %for cells
    n_r_ep = (n(x_int_p) - n(x_int_p-1))./(n(x_int_p+1) - n(x_int_p));    
    %for r_w, start at 2 to avoid sampling off the grid ... assign
    %value of -1 for first values (would be true for 0 Neumann BC)
    n_r_wp = (n(x_int_p(2:end)-1) - n(x_int_p(2:end)-2))./(n(x_int_p(2:end)) - n(x_int_p(2:end)-1));
    n_r_wp = [-1;n_r_wp];
    
    
    %for collagen
    rho_r_ep = (rho(x_int_p) - rho(x_int_p-1))./(rho(x_int_p+1) - rho(x_int_p));    
    %for r_w, start at 2 to avoid sampling off the grid ... assign
    %value of -1 for first values (would be true for 0 Neumann BC)
    rho_r_wp = (rho(x_int_p(2:end)-1) - rho(x_int_p(2:end)-2))./(rho(x_int_p(2:end)) - rho(x_int_p(2:end)-1));
    rho_r_wp = [-1;rho_r_wp];
    
    
    %replace nan, inf vals
    n_r_ep(isnan(n_r_ep)) = 1;
    n_r_wp(isnan(n_r_wp)) = 1;
    n_r_ep(isinf(n_r_ep)) = 100;
    n_r_wp(isinf(n_r_wp)) = 100;
    rho_r_ep(isnan(rho_r_ep)) = 1;
    rho_r_wp(isnan(rho_r_wp)) = 1;
    rho_r_ep(isinf(rho_r_ep)) = 100;
    rho_r_wp(isinf(rho_r_wp)) = 100;
    
    %compute A_posx
    Axp_n = A_pos(sigma(n_r_ep),sigma(n_r_wp),ve(x_int_p),vw(x_int_p),x_int_p,1);
    Axp_rho = A_pos(sigma(rho_r_ep),sigma(rho_r_wp),ve(x_int_p),vw(x_int_p),x_int_p,1);
    
    
    %Construct Aneg (negative matrix)
    
    %sensors for cells
    n_r_en = (n(x_int_n(1:end-1)+1) - n(x_int_n(1:end-1)+2))./(n(x_int_n(1:end-1)) - n(x_int_n(1:end-1)+1));
    n_r_en = [n_r_en;-1];
    n_r_wn = (n(x_int_n) - n(x_int_n+1))./(n(x_int_n-1) - n(x_int_n));

    %sensors for collagen
    rho_r_en = (rho(x_int_n(1:end-1)+1) - rho(x_int_n(1:end-1)+2))./(rho(x_int_n(1:end-1)) - rho(x_int_n(1:end-1)+1));
    rho_r_en = [rho_r_en;-1];
    rho_r_wn = (rho(x_int_n) - rho(x_int_n+1))./(rho(x_int_n-1) - rho(x_int_n));

    %replace nan, inf values
    n_r_en(isnan(n_r_en)) = 1;
    n_r_wn(isnan(n_r_wn)) = 1;
    n_r_en(isinf(n_r_en)) = 100;
    n_r_wn(isinf(n_r_wn)) = 100;
    rho_r_en(isnan(rho_r_en)) = 1;
    rho_r_wn(isnan(rho_r_wn)) = 1;
    rho_r_en(isinf(rho_r_en)) = 100;
    rho_r_wn(isinf(rho_r_wn)) = 100;
    
    
    %compute A_posx
    Axn_n = A_neg(sigma(n_r_en),sigma(n_r_wn),ve(x_int_n),vw(x_int_n),x_int_n,1);
    Axn_rho = A_neg(sigma(rho_r_en),sigma(rho_r_wn),ve(x_int_n),vw(x_int_n),x_int_n,1);
    
    %final derivative calculations
    dndt = (q(5)*D2 - 1/dx*(Axp_n + Axn_n))*n;
    
    drhodt = q(4)*n - 1/dx*(Axp_rho + Axn_rho)*rho;
    
    yprime = [dndt ; drhodt ; v];



end