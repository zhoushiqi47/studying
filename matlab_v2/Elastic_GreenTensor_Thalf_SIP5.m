function G=Elastic_GreenTensor_Thalf_SIP5(omega,kp,ks,y,x)
%% x2=0
%% Traction free wave Green Tensor
% the integral path is SIP
% 向量形式积分
%y is point source
% 2016 12 17 
%modify 2017 11 28
%matlab 默认branch Re>0

%rayleigh root precise expresion

[r1,r2]=raylei_eq(kp^2/(ks^2));
kr=ks*r1;
kr2=ks^2*r2;

n=size(x);
G=zeros(4,n(2));
L=5*kr;
mu = omega.^2./(ks.^2);
a=abs(x(1,:)-y(1,:));
ss=sign(x(1,:)-y(1,:));
b=y(2,:);
%compute rayleigh wave
mus=(ks^2-kr2)^(1/2);
mup=(kp^2-kr2)^(1/2);
beta=ks^2-2*kr2;
dedelta=8*kr*mus*mup-4*beta*kr-4*kr2*kr*mus/mup-4*kr2*kr*mup/mus;
Es=exp(1i*mus*b);
Ep=exp(1i*mup*b);

Gr1= 1i*(mus*beta*Es+2*mus*kr2*Ep)/(2*pi*dedelta*mu);
Gr2= 1i*(2*kr*mus*mup*Es-kr*beta*Ep)/(2*pi*dedelta*mu);
Gr3= 1i*(kr*beta*Es-2*kr*mus*mup*Ep)/(2*pi*dedelta*mu);
Gr4= 1i*(2*kr2*mup*Es+mup*beta*Ep)/(2*pi*dedelta*mu);

GR1= 1i*pi*Gr1.*exp(1i*kr*(a));
GR2= 1i*pi*Gr2.*exp(1i*kr*(a));
GR3= 1i*pi*Gr3.*exp(1i*kr*(a));
GR4= 1i*pi*Gr4.*exp(1i*kr*(a));

G1=integral(@(t) A1(t,a,b,kp,ks,mu,Gr1) ,-L,L, 'ArrayValued',true);
G2=integral(@(t) A2(t,a,b,kp,ks,mu,Gr2) ,-L,L, 'ArrayValued',true);
G3=integral(@(t) A3(t,a,b,kp,ks,mu,Gr3) ,-L,L, 'ArrayValued',true);
G4=integral(@(t) A4(t,a,b,kp,ks,mu,Gr4) ,-L,L, 'ArrayValued',true);



G(1,:)=G1+GR1;
G(2,:)=(G2+GR2).*ss;
G(3,:)=(G3+GR3).*ss;
G(4,:)=G4+GR4;

end


function f=A1(t,a,b,kp,ks,mu,c)
[r1,r2]=raylei_eq(kp^2/(ks^2));
kr=ks*r1;
kr2=ks^2*r2;
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
B=mus*beta;
C=2*mus*t^2;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=1i*(B*Es+C*Ep)/(2*pi*delta*mu)-(1/(t-kr)-1/(t+kr))*c;
end

function f=A2(t,a,b,kp,ks,mu,c)
[r1,r2]=raylei_eq(kp^2/(ks^2));
kr=ks*r1;
kr2=ks^2*r2;
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
B=2*t*mus*mup;
C=-t*beta;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=1i*(B*Es+C*Ep)/(2*pi*delta*mu)-(1/(t-kr)+1/(t+kr))*c;
end

function f=A3(t,a,b,kp,ks,mu,c)
[r1,r2]=raylei_eq(kp^2/(ks^2));
kr=ks*r1;
kr2=ks^2*r2;
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
C=-2*t*mus*mup;
B=t*beta;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=1i*(B*Es+C*Ep)/(2*pi*delta*mu)-(1/(t-kr)+1/(t+kr))*c;
end

function f=A4(t,a,b,kp,ks,mu,c)
[r1,r2]=raylei_eq(kp^2/(ks^2));
kr=ks*r1;
kr2=ks^2*r2;
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
B=2*t^2*mup;
C=mup*beta;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=1i*(B*Es+C*Ep)/(2*pi*delta*mu)-(1/(t-kr)-1/(t+kr))*c;
end