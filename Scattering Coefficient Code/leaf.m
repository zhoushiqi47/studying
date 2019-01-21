function [x,y]=leaf(t,type)
%% p-leaf
 theta = t;
 rho = (1+cos(theta));
 drho = -sin(theta);
 ddrho =-cos(theta);
 if type==1
     x = rho.*cos(theta);   %% primal function
     y = rho.*sin(t);
 else if type==2
         x = -rho.*sin(t)+drho.*cos(t);  %% derivative of order one
         y = rho.*cos(t)+sin(t).*drho ;
     else if type==3     %% derivative of order two
             x = ddrho.*cos(t) - 2*drho.*sin(t) + rho.*cos(t) ;
             y = ddrho.*sin(t) + 2*drho.*cos(t) - rho.*sin(t);
         end
     end
     
 end
 return