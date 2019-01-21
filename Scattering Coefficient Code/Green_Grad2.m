function U=Green_Grad2(k,x1,x2,dx1,dx2s)
%% The calculation of the double gradient
z1=z(1,:)-z0(1,:);
z2=z(2,:)-z0(2,:);
d=sqrt(z1.*z1+z2.*z2);
vad1= besselh(0,k*d);
vad2= besselh(2,k*d);
U=z1(1,:)*0;
U = 1i/4*k*[val.*z1./d; val.*z2./d]
 

