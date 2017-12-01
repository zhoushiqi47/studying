function G=TractionDGreenTensor_2D2(omega,kp,ks,y,x)
%% n=(0,1)
%% operator T to x
%% compute the traction green tensor of elastic wave operator in 2D
%% shear wave velocity
%% x2=0
%% pressure wave Green Tensor
%% modify 2017 11 30



% 向量形式积分
%y is point source
% 2016 12 16 
tic
n=size(x);
G=zeros(4,n(2));

a=x(1,:)-y(1,:);
b=y(2,:);

l=5*ks;
k1=kp;
k2=ks;

% 设置间接变量

Gc1=integral(@(t) A1(t,a,b,k1,k2) ,-l,l, 'ArrayValued',true);
Gc2=integral(@(t) A2(t,a,b,k1,k2) ,-l,l, 'ArrayValued',true);
Gc3=integral(@(t) A3(t,a,b,k1,k2) ,-l,l, 'ArrayValued',true);
Gc4=integral(@(t) A4(t,a,b,k1,k2) ,-l,l, 'ArrayValued',true);

G(1,:)=Gc1;
G(2,:)=Gc2;
G(3,:)=Gc3;
G(4,:)=Gc4;
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


