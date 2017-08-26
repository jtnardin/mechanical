function yprime = mechanical_ode_TM_base_2d(t,y,q,dx,dy,D1X,D1Y,D2X,D2Y,D2XY,...
    D2_bd_D,D2Xbd_N,D2Ybd_N,A_pos,A_pos_0,A_pos_1,A_neg,A_neg_0,A_neg_1,...
    xy_int,bd,x_bd_0,x_bd_l,y_bd_0,y_bd_l,total,x_int,y_int,xn,yn)

%       t


    x = linspace(-1,1,xn);
    y1 = linspace(0,1,yn);
    [Y,X] = meshgrid(y1,x);



    n = y(1:total);
    rho = y(total+1:2*total);
    u = y(2*total+1:3*total);
    
    
    
    v = velocity_solve(q,D1X,D1Y,D2X,D2Y,D2XY,D2_bd_D,n,rho,repmat(u,2,1),xy_int,bd,total);
    
    
    vx = v(1:total);
    vy = v(total+1:end);
    
    %%%% don't forget boundary conditions (8/25)
    
    ve_int = vx(x_int+1);
    ve_x_0 = vx(x_bd_0+1);
    ve_x_1 = vx(x_bd_l);
    % do we need a east velocity at x = 1?
        
    vw_int = vx(x_int-1);
    vw_x_1 = vx(x_bd_l-1);
    vw_x_0 = vx(x_bd_0);
    % do we need a west velocity at x = 1?
    
    vn_int = vy(y_int+xn);
    vn_y_0 = vy(y_bd_0+xn);
    vn_y_1 = vy(y_bd_l);
    % do we need a north velocity at y = 1?
        
    vs_int = vy(y_int-xn);
    vs_y_1 = vy(y_bd_l-xn);
    vs_y_0 = vy(y_bd_0);
    % do we need a south velocity at y = 0?
        
    %estimate positive, negative velocity locations
    Vx_pos_loc = ve_int + vw_int >= 0;
    Vx_neg_loc = ve_int + vw_int <  0;

    Vy_pos_loc = vn_int + vs_int >= 0;
    Vy_neg_loc = vn_int + vs_int < 0;

    %compute velocities at neg, pos points
    vw_pos = 1/dx*vw_int(Vx_pos_loc);
    ve_pos = 1/dx*ve_int(Vx_pos_loc);

    vw_neg = 1/dx*vw_int(Vx_neg_loc);
    ve_neg = 1/dx*ve_int(Vx_neg_loc);

    vs_pos = 1/dy*vs_int(Vy_pos_loc);
    vn_pos = 1/dy*vn_int(Vy_pos_loc);

    vs_neg = 1/dy*vs_int(Vy_neg_loc);
    vn_neg = 1/dy*vn_int(Vy_neg_loc);


    %sigma for flux limiters
    sigma = @(r) (r+abs(r))./(1+abs(r));

    
    
    %Construction of Ax involves splitting into positive, negative components    
    
    %locate x_int points for positive x vel
    x_intp = x_int(Vx_pos_loc); 
    %which points are on the second row?
    x_intp_row1 = mod(x_intp,xn)==2;
    %locate x_int points for negative x vel
    x_intn = x_int(Vx_neg_loc);
    %which points are in penultimate row?
    x_intn_rown = mod(x_intn,xn)==(xn-1);
   
    [n_ep,n_wp] = positive_sensor(n,x_intp,x_intp_row1,1);
    [rho_ep,rho_wp] = positive_sensor(rho,x_intp,x_intp_row1,1);
    [n_en,n_wn] = negative_sensor(n,x_intn,x_intn_rown,1);
    [rho_en,rho_wn] = negative_sensor(rho,x_intn,x_intn_rown,1);
    
    
    
    %compute A_posx
    Axp_n = A_pos(sigma(n_ep),sigma(n_wp),ve_pos,vw_pos,x_intp,1);
    Axp_rho = A_pos(sigma(rho_ep),sigma(rho_wp),ve_pos,vw_pos,x_intp,1);
        
    %Construct Aneg (negative matrix)
       
    %compute A_negx
    Axn_n = A_neg(sigma(n_en),sigma(n_wn),ve_neg,vw_neg,x_intn,1);
    Axn_rho = A_neg(sigma(rho_en),sigma(rho_wn),ve_neg,vw_neg,x_intn,1);
    
        
    %%% construction of Ay
    
    
    %locate y_int points for positive y vel
    y_intp = y_int(Vy_pos_loc); 
    %which points are on the second column?
    y_intp_col1 = y_intp <= 2*xn;
    %locate xy_int points for negativ x vel
    y_intn = y_int(Vy_neg_loc);
    %how many interior points are in penultimate column
    y_intn_coln = y_intn >= (yn-2)*xn;
   
    [n_np,n_sp] = positive_sensor(n,y_intp,y_intp_col1,xn);
    [rho_np,rho_sp] = positive_sensor(rho,y_intp,y_intp_col1,xn);
    [n_nn,n_sn] = negative_sensor(n,y_intn,y_intn_coln,xn);
    [rho_nn,rho_sn] = negative_sensor(rho,y_intn,y_intn_coln,xn);
  
    
    
    %compute A_posy
    Ayp_n = A_pos(sigma(n_np),sigma(n_sp),vn_pos,vs_pos,y_intp,xn);
    Ayp_rho = A_pos(sigma(rho_np),sigma(rho_sp),vn_pos,vs_pos,y_intp,xn);
    
    
    %Construct Aneg (negative matrix)
    %compute A_posx
    Ayn_n = A_pos(sigma(n_nn),sigma(n_sn),vn_neg,vs_neg,y_intn,xn);
    Ayn_rho = A_pos(sigma(rho_nn),sigma(rho_sn),vn_neg,vs_neg,y_intn,xn);
    
    
    
    %final derivative calculations
    dndt = ( q(9)*(D2X(xy_int,total) + D2Y(xy_int,total)  + D2Xbd_N(x_bd_0,x_bd_l,total)...
        + D2Ybd_N(y_bd_0,y_bd_l,total)) - (Axp_n + Axn_n + Ayp_n + Ayn_n) + q(8)*sparse(1:total,1:total,1-n,total,total))*n;
    
    drhodt =  -(Axp_rho + Axn_rho + Ayp_rho + Ayn_rho )*rho;
    
    yprime = [dndt ; drhodt ; v];



end