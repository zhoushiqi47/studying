function [x,y]=Penut(t,type)
r = 1;
 if type==1
   x = r*(cos(t) + 0.2*cos(3*t));
    y = r*(sin(t) + 0.2*sin(3*t)) +10;
else if type==2
        x = -r*(sin(t) + 3*0.2*sin(3*t));
        y = r*(cos(t) + 3*0.2*cos(3*t));
   else if type==3
            x = - r*( cos(t) + 3*3*0.2*cos(3*t));
             y = - r* ( sin(t) + 3*3*0.2*sin(3*t));
       else if type ==4
               x =  r*( sin(t) + 3*3*3*0.2*sin(3*t));
               y = - r* ( cos(t) + 3*3*3*0.2*sin(3*t));
           end
        end
    end
 end
 return 