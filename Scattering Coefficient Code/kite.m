function [x,y]=kite(t,type)

if type==1
    x = cos(t)+0.65*cos(2*t)-0.65 ;   %% primal function
    y = 1.5*sin(t) ;
else if type==2
        x = -sin(t)-1.3*sin(2*t);  %% derivative of order one
        y = 1.5*cos(t);
    else if type==3     %% derivative of order two
            x = -cos(t)-2.6*cos(2*t);
            y = -1.5*sin(t);
        else if type==4
            x = sin(t) + 5.2*sin(2*t);
            y = -1.5*cos(t);
            end
        end
    end
end
% x = 2*x;
% y = 2*y;

return 
            