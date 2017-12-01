function G=TractionDGreenTensor_2D_prin(omega,kp,ks,y,x)
%% n=(0,1)
%% operator T to x
%% compute the traction green tensor of elastic wave operator in 2D
%% shear wave velocity
%% x2=0
%% pressure wave Green Tensor
%% 2017 11 30
%% principle term



% 向量形式积分
%y is point source
% 2016 12 16 
n=size(x);
G=zeros(4,n(2));
c=sign(x(1,:)-y(1,:));
a=abs(x(1,:)-y(1,:));
b=y(2,:);
d=sqrt(a.^2+b.^2);

ts=ks*(a)./d;
tp=kp*(a)./d;

L2=3*ks;
L1=ks;

smus=(ks^2-ts.^2).^(1/2);
smup=(kp^2-ts.^2).^(1/2);
sdelta=ts.^2+smus.*smup;

pmus=(ks^2-tp.^2).^(1/2);
pmup=(kp^2-tp.^2).^(1/2);
pdelta=tp.^2+pmus.*pmup;


Es=ks/sqrt(2*pi)./sdelta.*b./d./sqrt(ks*d).*exp(1i*ks*d-1i*pi/4);
Ep=kp/sqrt(2*pi)./pdelta.*b./d./sqrt(kp*d).*exp(1i*kp*d-1i*pi/4);

Gp1=smus.*smup.*Es+tp.^2.*Ep;
Gp2=(ts.*smus.*Es-tp.*pmus.*Ep).*c;
Gp3=(ts.*smup.*Es-tp.*pmup.*Ep).*c;
Gp4=ts.^2.*Es+pmus.*pmup.*Ep;

%Gc1=integral(@(t) A1(t,x(:,:),y(:,:),kp,ks) ,L1,L2, 'ArrayValued',true)+integral(@(t) A1(t,x(:,:),y(:,:),kp,ks) ,-L2,-L1, 'ArrayValued',true);
%Gc2=integral(@(t) A2(t,x(:,:),y(:,:),kp,ks) ,L1,L2, 'ArrayValued',true)+integral(@(t) A2(t,x(:,:),y(:,:),kp,ks) ,-L2,-L1, 'ArrayValued',true);
%Gc3=integral(@(t) A3(t,x(:,:),y(:,:),kp,ks) ,L1,L2, 'ArrayValued',true)+integral(@(t) A3(t,x(:,:),y(:,:),kp,ks) ,-L2,-L1, 'ArrayValued',true);
%Gc4=integral(@(t) A4(t,x(:,:),y(:,:),kp,ks) ,L1,L2, 'ArrayValued',true)+integral(@(t) A4(t,x(:,:),y(:,:),kp,ks) ,-L2,-L1, 'ArrayValued',true);

G(1,:)=Gp1;
G(2,:)=Gp2;
G(3,:)=Gp3;
G(4,:)=Gp4;

end

function f=A1(t,x,y,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
delta=t^2+mus*mup;
B=mus*mup;
C=t*t;
Es=exp(1i*mus*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
Ep=exp(1i*mup*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
f=(B*Es+C*Ep)/(2*pi*delta);
end

function f=A2(t,x,y,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
delta=t^2+mus*mup;
B=t*mus;
Es=exp(1i*mus*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
Ep=exp(1i*mup*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
f=(B*Es-B*Ep)/(2*pi*delta);
end

function f=A3(t,x,y,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
delta=t^2+mus*mup;
B=t*mup;
Es=exp(1i*mus*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
Ep=exp(1i*mup*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
f=(B*Es-B*Ep)/(2*pi*delta);
end

function f=A4(t,x,y,kp,ks)
mus=(ks^2-t^2)^(1/2);
mup=(kp^2-t^2)^(1/2);
delta=t^2+mus*mup;
B=mus*mup;
C=t*t;
Es=exp(1i*mus*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
Ep=exp(1i*mup*(y(2,:))+1i*t*(x(1,:)-y(1,:)));
f=(C*Es+B*Ep)/(2*pi*delta);
end


