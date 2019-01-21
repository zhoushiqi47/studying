function [x,y]=p_leaf(t,type) 
%rotate
 a=sin(pi/4);
 p = 4;
 d = 1;
 rho = d*(1+0.2*cos(p*t));
 drho = -d*p*0.2*sin(p*t);
 ddrho =-d*p^2*0.2*cos(p*t);
 if type==1
     x1 = rho.*cos(t);   %% primal function
     y1 = rho.*sin(t);
     x = a*x1-a*y1;
     y = a*x1+a*y1+10;
 else if type==2
         x1 = -rho.*sin(t)+drho.*cos(t);  %% derivative of order one
         y1 = rho.*cos(t)+sin(t).*drho ;
         x = a*x1-a*y1;
         y = a*x1+a*y1;
     else if type==3     %% derivative of order two
             x1 = ddrho.*cos(t) - 2*drho.*sin(t) + rho.*cos(t) ;
             y1 = ddrho.*sin(t) + 2*drho.*cos(t) - rho.*sin(t);
             x = a*x1-a*y1;
             y = a*x1+a*y1;
         end
     end
 end
 return