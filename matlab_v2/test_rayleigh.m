omega=2*pi;
ks=4*pi;
kp=2*pi;
tau=kp^2/(ks^2);

[r1,r2]=raylei_eq(tau);
kr=ks*r1;
kr2=ks^2*r2;



mu = omega.^2./(ks.^2);
a=0.1;

b=0.1;

musr=(ks^2-kr2)^(1/2);
mupr=(kp^2-kr2)^(1/2);
betar=ks^2-2*kr2;
dedeltar=8*kr*musr*mupr-4*betar*kr-4*kr2*kr*musr/mupr-4*kr2*kr*mupr/musr;
Esr=exp(1i*musr*b);
Epr=exp(1i*mupr*b);

Gr1= 1i*(musr*betar*Esr+2*musr*kr2*Epr)/(2*pi*dedeltar*mu);

t=linspace(13.47,13.48,5000);
mus=(ks^2-t.^2).^(1/2);
mup=(kp^2-t.^2).^(1/2);
beta=ks^2-2*t.^2;
delta=beta.^2+4*t.^2.*mus.*mup;
ndelta=1./delta;
B=mus.*beta;
C=2*mus.*t.^2;
Es=exp(1i*mus*b+1i*t*a);
Ep=exp(1i*mup*b+1i*t*a);

f1=1i*(B.*Es+C.*Ep)./(2*pi*delta*mu);
f2=(1./(t-kr)-1./(t+kr))*1i*(musr*betar*Esr+2*musr*kr2*Epr)/(2*pi*dedeltar*mu);

f=f1-f2;

figure;
plot(t,real(1.0./delta))


