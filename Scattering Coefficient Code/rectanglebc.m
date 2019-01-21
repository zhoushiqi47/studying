function [x,y]=rectanglebc(t,type)
%% the rectangle shape %%%%%%%%%%%%%%%%%%
 radius = 0.5;
if type==1
    x = radius*( cos(t).^3 + cos(t) );   %% primal function
    y = radius*( sin(t).^3 + sin(t) ) ;
else if type==2
        x = -radius*sin(t).*(  3*cos(t).^2 + 1  );  %% derivative of order one
        y = radius*cos(t).*(   3*sin(t).^2 + 1 );
    else if type==3     %% derivative of order two
            x = -radius*( 3*cos(t).^3 + cos(t) -6*cos(t).*sin(t).^2 );
            y =  radius*( -(3*sin(t).^3 + sin(t)) + 6*sin(t).*cos(t).^2 ) ;
        end
    end
end

return
