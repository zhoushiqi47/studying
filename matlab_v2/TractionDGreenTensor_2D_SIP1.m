function G=TractionDGreenTensor_2D_SIP1(omega,kp,ks,y,x)
%% n=(0,1)
%% operator T to x
%% compute the traction green tensor of elastic wave operator in 2D
%% shear wave velocity
%% x2=0
%% pressure wave Green Tensor



%The integral path is SIP
% 向量形式积分
%y is point source
% 2016 12 16 
%modify 2017 11 28 new integral path
tic
n=size(x);
G=zeros(4,n(2));

L=10*(kp+ks);

signs=sign(x(1,:)-y(1,:));
a=abs(x(1,:)-y(1,:));
b=y(2,:);
% 设置间接变量


G1=integral(@(t) A1(t,a,b,kp,ks) ,-ks,ks, 'ArrayValued',true);
G2=integral(@(t) A2(t,a,b,kp,ks) ,-ks,ks, 'ArrayValued',true);
G3=integral(@(t) A3(t,a,b,kp,ks) ,-ks,ks, 'ArrayValued',true);
G4=integral(@(t) A4(t,a,b,kp,ks) ,-ks,ks, 'ArrayValued',true);

Gc1=integral(@(t) A1(t,a,b,kp,ks) ,-L+L/2*1i,-ks+0i, 'ArrayValued',true);
Gc2=integral(@(t) A2(t,a,b,kp,ks) ,-L+L/2*1i,-ks+0i, 'ArrayValued',true);
Gc3=integral(@(t) A3(t,a,b,kp,ks) ,-L+L/2*1i,-ks+0i, 'ArrayValued',true);
Gc4=integral(@(t) A4(t,a,b,kp,ks) ,-L+L/2*1i,-ks+0i, 'ArrayValued',true);

G(1,:)=G1+2*real(Gc1);
G(2,:)=(G2+2*real(Gc2)).*signs;
G(3,:)=(G3+2*real(Gc3)).*signs;
G(4,:)=G4+2*real(Gc4);
toc
end

function f=A1(t,a,b,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
delta=t^2+mus*mup;
B=mus*mup;
C=t*t;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=(B*Es+C*Ep)/(2*pi*delta);
end

function f=A2(t,a,b,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
delta=t^2+mus*mup;
B=t*mus;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=(B*Es-B*Ep)/(2*pi*delta);
end

function f=A3(t,a,b,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
delta=t^2+mus*mup;
B=t*mup;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=(B*Es-B*Ep)/(2*pi*delta);
end

function f=A4(t,a,b,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
delta=t^2+mus*mup;
B=mus*mup;
C=t*t;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);
f=(C*Es+B*Ep)/(2*pi*delta);
end


