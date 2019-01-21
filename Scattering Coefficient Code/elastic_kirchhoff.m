clear;

lamda=1/2;

mu=1/4;
omega=4*pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);
kappa = kp/ks;
lambda= lamda;
kappa1 = ks/kp;
n = 256;
node = 0:2*n-1;
t = pi*node(:)/n;

bctype = 1;

if bctype==1
    [x1,x2]=circlebc(t,1);   
    [dx1,dx2]=circlebc(t,2);
    [ddx1,ddx2]=circlebc(t,3);

else if bctype==2
    [x1,x2]=Penut(t,1);   
    [dx1,dx2]=Penut(t,2);
    [ddx1,ddx2]=Penut(t,3);
    else if bctype==3
            [x1,x2]=p_leaf(t,1);   
            [dx1,dx2]=p_leaf(t,2);
            [ddx1,ddx2]=p_leaf(t,3);
        else
            [x1,x2]=myrectangle(t,1);   
            [dx1,dx2]=myrectangle(t,2);
            [ddx1,ddx2]=myrectangle(t,3);
        end
    end
end
distance = sqrt( dx1.*dx1+dx2.*dx2 );
trp1 = zeros(2*n,1);
trp2 = zeros(2*n,1);
trs1 = zeros(2*n,1);
trs2 = zeros(2*n,1);

t_x1 = dx1./distance;
t_x2 = dx2./distance;
t_x = [t_x1, t_x2];
n_x1 = t_x2;
n_x2 = -t_x1;
n_x =[n_x1, n_x2];
R1=zeros(2,2*n);
R2=zeros(2,2*n);

theta=0;
p0 = [sin(theta);cos(theta)];
pn0 = [cos(theta);-sin(theta)];
%% compute traction of total field in plane
for j=1:2*n
    nj = n_x(j,:);
    alpha0=nj*p0;
    
    if alpha0>=0
        R1(:,j)=0;
    else
    tj = t_x(j,:);
    beta = kappa*alpha0 - sqrt(kappa^2*alpha0^2-(kappa^2-1));
    p1 = p0 -2*alpha0*(nj.');
    p2 = kappa*p0- beta*(nj.');
    pn2 = [p2(2);-p2(1)];
    A0 = 2*kappa*alpha0^2-kappa-alpha0*beta;
    A1 = kappa - alpha0*beta;
    A2 = 2*alpha0*(tj*p0);
    T0 = 1i*kp*(lambda*nj.'+2*mu*alpha0*p0);
    T1 = 1i*kp*(lambda*nj.'+2*mu*(nj*p1)*p1)*A1/A0;
    T2 = 1i*ks*mu*((nj*p2)*pn2+(nj*pn2)*p2)*A2/A0;
    R1(:,j) = (T0+T1+T2)/kp;
    
    gamma = kappa1*alpha0 - sqrt(kappa1^2*alpha0^2-(kappa1^2-1));
    q0 = p0;
    qn0 = pn0;
    q1 = kappa1*q0 - gamma*(nj.');
    q2 = q0 - 2*alpha0*(nj.');
    qn2 = [q2(2); -q2(1)];
    A0 = 2*kappa1*alpha0^2-kappa1-alpha0*gamma;
    A1 = -2*alpha0*(tj*q0);
    A2 = kappa1- gamma*alpha0;
    T0 = 1i*ks*mu*(alpha0*qn0+(nj*qn0)*q0);
    T1 = 1i*kp*(lambda*nj.'+2*mu*(nj*q1)*q1)*A1/A0;
    T2 = 1i*ks*mu*((nj*q2)*qn2+(nj*qn2)*q2)*A2/A0;
    R2(:,j) = (T0+T1+T2)/ks;
    end
end



 A = singlelayar_fullspace(omega, lambda, mu, n, bctype);
 
 
 
 fp1 = sin(theta)*exp(1i*kp*(x1*sin(theta)+x2*cos(theta)));
 fp2 = cos(theta)* exp(1i*kp*(x1*sin(theta)+x2*cos(theta)));
 fp = [fp1;fp2];
 fs1 =cos(theta)* exp(1i*ks*(x1*sin(theta)+x2*cos(theta)));
 fs2 =- sin(theta)*exp(1i*ks*(x1*sin(theta)+x2*cos(theta)));
 fs = [fs1;fs2];

 
 f = [fp,fs];
 
 Phi = A\f;
 
 R= Phi;
 x=linspace(0,2,2*n);
 figure; 
  subplot(2,2,1);plot(x,imag(R(1:2*n,1)/kp./exp(1i*kp*(x1*sin(theta)+x2*cos(theta)))));
  hold on ; plot(x,imag(R1(1,:)),'r');
  subplot(2,2,2) ;plot(x,imag(R1(1,:)));
  
  subplot(2,2,3);plot(x,imag(R(2*n+1:4*n,1)/kp./exp(1i*kp*(x1*sin(theta)+x2*cos(theta)))));
  hold on ; plot(x,imag(R1(2,:)),'r');
  subplot(2,2,4); plot(x,imag(R1(2,:)));
  
   figure; 
  subplot(2,2,1);plot(x,imag(R(1:2*n,2)/ks./exp(1i*ks*(x1*sin(theta)+x2*cos(theta)))));
  hold on ; plot(x,imag(R2(1,:)),'r');
  subplot(2,2,2) ;plot(x,imag(R2(1,:)));
  
  subplot(2,2,3);plot(x,imag(R(2*n+1:4*n,2)/ks./exp(1i*ks*(x1*sin(theta)+x2*cos(theta)))));
  hold on ; plot(x,imag(R2(2,:)),'r');
  subplot(2,2,4); plot(x,imag(R2(2,:)));


  



  
 