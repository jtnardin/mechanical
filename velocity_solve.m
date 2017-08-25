%function velocity_solve.m written 8-24-17 by JTN to solve a 2d velocity
%field based on Tranquillo Murray (1992)

%8-24-17 : can solve for interior velocity ; be sure to incorprate boundary
%conditions soon

function v = velocity_solve(q,D1X,D1Y,D2X,D2Y,D2XY,D2_bd,n,rho,u,xy_int,bd,total)


      xy_inty = xy_int + total;
      bdy = bd + total;

      LHS = (q(1)+q(2))*D2X(xy_int,2*total) + q(2)/2*D2Y(xy_int,2*total) + q(1)/2*D2X(xy_inty,2*total) +...
          (q(1)+q(2))*D2Y(xy_inty,2*total) + (q(1)/2+q(2))*D2XY(xy_int,xy_inty,2*total) + ...
          (q(1)/2+q(2))*D2XY(xy_inty,xy_int,2*total) + D2_bd(bd,2*total) + D2_bd(bdy,2*total);

     
      RHS = q(5)*repmat(rho,2,1).*u - ((q(3)/(1-q(4))+q(4)/(1-2*q(4)))*D2X(xy_int,2*total) + ...
          q(3)/2/(1-q(4))*D2Y(xy_int,2*total) + q(3)/2/(1-q(4))*D2X(xy_inty,2*total) +...
          (q(3)/(1-q(4))+q(4)/(1-2*q(4)))*D2Y(xy_inty,2*total) + (q(3)/2/(1-q(4)) +...
          q(4)/(1-2*q(4)))*D2XY(xy_int,xy_inty,2*total) + (q(3)/2/(1-q(4)) + ...
          q(4)/(1-2*q(4)))*D2XY(xy_inty,xy_int,2*total))*u - q(6)*(D1X(xy_int,2*total) +...
          D1Y(xy_inty,2*total))*repmat(rho.*n./(1+q(7)*n.^2),2,1);
      
%     v = -1/q(1)*(D2 + D2bd)\(q(2)*((D1*rho).*(D1*u) + (rho).*(D2*u)) + ...
%         q(3)*((D1*rho).*n + rho.*(D1*n)) - q(6)*rho.*n);


      v = LHS\RHS;







end