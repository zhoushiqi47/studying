function [x,y]=kite(t,type)
r = 1;
 if type==1
   x = r*(cos(t) + 0.65*cos(2*t))-0.65;
    y = r*1.5*sin(t) +10;
else if type==2
        x = -r*(sin(t) + 2*0.65*sin(2*t));
        y = r*1.5*cos(t);
   else if type==3
            x = - r*( cos(t) + 2* 2*0.65*cos(2*t));
             y = - r* 1.5* sin(t) ;
       else if type==3
               x=  r*( sin(t) + 2*2* 2*0.65*sin(2*t));
               y = - r* 1.5*cos(t) ;
           end
        end
    end
 end
 return 