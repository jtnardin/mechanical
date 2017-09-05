%function velocity_solve.m written 8-24-17 by JTN to solve a 2d velocity
%field based on Tranquillo Murray (1992)

%8-24-17 : can solve for interior velocity ; be sure to incorprate boundary
%conditions soon

function v = velocity_solve(q,D1X,D1Y,D2X,D2X0,D2XL,D2Y,D2Y0,D2YL,D2XY,n,...
    rho,u,xy_int,x_int,y_int,x_bd_0,x_bd_l,y_bd_0,y_bd_l,total)

        %need D2X0, D2XL, D2Y0, D2YL, x_int , x_bd_0, x_bd_l , y_bd_0 ,
        %y_bd_l


      xy_inty  = xy_int + total;
      x_inty   = x_int + total;
      y_inty   = y_int + total;
      x_bd_0y  = x_bd_0 + total;
      x_bd_ly  = x_bd_l + total;
      y_bd_0y  = y_bd_0 + total;
      y_bd_ly  = y_bd_l + total;
      

      LHS = (q(1)+q(2))*(D2X(x_int,2*total) + D2X0(x_bd_0,2*total) +...
          D2XL(x_bd_l,2*total)) + q(1)/2*(D2Y(y_int,2*total) + ...
          D2Y0(y_bd_0,2*total) + D2YL(y_bd_l,2*total)) + q(1)/2*(D2X(x_inty,2*total) + ...
          D2X0(x_bd_0y,2*total) + D2XL(x_bd_ly,2*total)) + (q(1)+q(2))*(D2Y(y_inty,2*total) + ...
          D2Y0(y_bd_0y,2*total) + D2YL(y_bd_ly,2*total)) + (q(1)/2+q(2))*D2XY(xy_int,xy_inty,2*total) + ...
          (q(1)/2+q(2))*D2XY(xy_inty,xy_int,2*total);

     
      RHS = q(5)*repmat(rho,2,1).*u - ((q(3)/(1-q(4))+q(4)/(1-2*q(4)))*(D2X(x_int,2*total) + ...
          D2X0(x_bd_0,2*total) + D2XL(x_bd_l,2*total)) + q(3)/2/(1-q(4))*(D2Y(y_int,2*total) + ...
          D2Y0(y_bd_0,2*total) + D2YL(y_bd_l,2*total)) + q(3)/2/(1-q(4))*(D2X(x_inty,2*total) + ...
          D2X0(x_bd_0y,2*total) + D2XL(x_bd_ly,2*total)) + (q(3)/(1-q(4))+q(4)/(1-2*q(4)))*(D2Y(y_inty,2*total) + ...
          D2Y0(y_bd_0y,2*total) + D2YL(y_bd_ly,2*total)) + (q(3)/2/(1-q(4)) +...
          q(4)/(1-2*q(4)))*D2XY(xy_int,xy_inty,2*total) + (q(3)/2/(1-q(4)) + ...
          q(4)/(1-2*q(4)))*D2XY(xy_inty,xy_int,2*total))*u - q(6)*(D1X(xy_int,2*total) +...
          D1Y(xy_inty,2*total))*repmat(rho.*n./(1+q(7)*n.^2),2,1);
      
%     v = -1/q(1)*(D2 + D2bd)\(q(2)*((D1*rho).*(D1*u) + (rho).*(D2*u)) + ...
%         q(3)*((D1*rho).*n + rho.*(D1*n)) - q(6)*rho.*n);


      v = LHS\RHS;







end