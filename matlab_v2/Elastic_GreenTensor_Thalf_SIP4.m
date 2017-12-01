function G=Elastic_GreenTensor_Thalf_SIP4(omega,kp,ks,y,x)
%% x2=0
%% Traction free wave Green Tensor
% the integral path is SIP
% 向量形式积分
%y is point source
% 2016 12 17 
%modify 2017 11 28
%matlab 默认branch Re>0
tau=kp/ks;
%rayleigh root precise expresion
tr=((16*tau^2 - 24)^2/(9*(16*tau^2 - 16)^2) + 8/(3*(16*tau^2 - 16)))/((16*tau^2 - 24)^3/(27*(16*tau^2 - 16)^3) + (((16*tau^2 - 24)^3/(27*(16*tau^2 - 16)^3) + (64*tau^2 - 96)/(3*(16*tau^2 - 16)^2) - 1/(32*tau^2 - 32))^2 - ((16*tau^2 - 24)^2/(9*(16*tau^2 - 16)^2) + 8/(48*tau^2 - 48))^3)^(1/2) + (64*tau^2 - 96)/(3*(16*tau^2 - 16)^2) - 1/(32*tau^2 - 32))^(1/3) + ((16*tau^2 - 24)^3/(27*(16*tau^2 - 16)^3) + (((16*tau^2 - 24)^3/(27*(16*tau^2 - 16)^3) + (4*(16*tau^2 - 24))/(3*(16*tau^2 - 16)^2) - 1/(2*(16*tau^2 - 16)))^2 - ((16*tau^2 - 24)^2/(9*(16*tau^2 - 16)^2) + 8/(3*(16*tau^2 - 16)))^3)^(1/2) + (4*(16*tau^2 - 24))/(3*(16*tau^2 - 16)^2) - 1/(2*(16*tau^2 - 16)))^(1/3) + (16*tau^2 - 24)/(3*(16*tau^2 - 16));
kr=sqrt(tr)*ks;
kr2=tr*ks^2;

n=size(x);
G=zeros(4,n(2));
L=10*(kp+ks);
mu = omega.^2./(ks.^2);

signs=sign(x(1,:)-y(1,:));
a=abs(x(1,:)-y(1,:));
b=y(2,:);
%compute rayleigh wave
mus=(ks^2-kr2)^(1/2);
mup=(kp^2-kr2)^(1/2);
beta=ks^2-2*kr2;
dedelta=8*kr*mus*mup-4*beta*kr-4*kr2*kr*mus/mup-4*kr2*kr*mup/mus;
Es=exp(1i*mus.*b+1i*kr.*a);
Ep=exp(1i*mup.*b+1i*kr.*a);

Gr1= -(mus*beta*Es+2*mus*kr2*Ep)/(dedelta*mu);
Gr2= -(2*kr*mus*mup*Es-kr*beta*Ep)/(dedelta*mu);
Gr3= -(kr*beta*Es-2*kr*mus*mup*Ep)/(dedelta*mu);
Gr4= -(2*kr2*mup*Es+mup*beta*Ep)/(dedelta*mu);

G1=integral(@(t) A1(t,a,b,kp,ks,mu) ,-ks,ks, 'ArrayValued',true);
G2=integral(@(t) A2(t,a,b,kp,ks,mu) ,-ks,ks, 'ArrayValued',true);
G3=integral(@(t) A3(t,a,b,kp,ks,mu) ,-ks,ks, 'ArrayValued',true);
G4=integral(@(t) A4(t,a,b,kp,ks,mu) ,-ks,ks, 'ArrayValued',true);

Gc1=integral(@(t) A1(t,a,b,kp,ks,mu) ,-L+L/2*1i,-ks+0i, 'ArrayValued',true);
Gc2=integral(@(t) A2(t,a,b,kp,ks,mu) ,-L+L/2*1i,-ks+0i, 'ArrayValued',true);
Gc3=integral(@(t) A3(t,a,b,kp,ks,mu) ,-L+L/2*1i,-ks+0i, 'ArrayValued',true);
Gc4=integral(@(t) A4(t,a,b,kp,ks,mu) ,-L+L/2*1i,-ks+0i, 'ArrayValued',true);

G(1,:)=G1+2*real(Gc1)+Gr1;
G(2,:)=(G2+2*real(Gc2)+Gr2).*signs;
G(3,:)=(G3+2*real(Gc3)+Gr3).*signs;
G(4,:)=G4+2*real(Gc4)+Gr4;

end


function f=A1(t,a,b,kp,ks,mu)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
B=mus*beta;
C=2*mus*t^2;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=1i*(B*Es+C*Ep)/(2*pi*delta*mu);
end

function f=A2(t,a,b,kp,ks,mu)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
B=2*t*mus*mup;
C=-t*beta;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=1i*(B*Es+C*Ep)/(2*pi*delta*mu);
end

function f=A3(t,a,b,kp,ks,mu)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
C=-2*t*mus*mup;
B=t*beta;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=1i*(B*Es+C*Ep)/(2*pi*delta*mu);
end

function f=A4(t,a,b,kp,ks,mu)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
beta=ks^2-2*t^2;
delta=beta^2+4*t^2*mus*mup;
B=2*t^2*mup;
C=mup*beta;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=1i*(B*Es+C*Ep)/(2*pi*delta*mu);
end