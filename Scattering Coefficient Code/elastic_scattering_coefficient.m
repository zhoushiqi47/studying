clear;

lamda=1/2;

mu=1/4;
f =8;
omega=2*f*pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);
kappa = kp/ks;
lambda= lamda;
kappa1 = ks/kp;

n1=128; %% discretization of incidence angle
theta = 2*pi/n1*(0:n1-1);
p0 = [cos(theta); sin(theta)]';
pn0 = [sin(theta); -cos(theta)]';


n2 = 512;  %% discretization of the boundary of obstacle
t = 2*pi/n2*(0:n2-1);

bctype = 2;

if bctype==1
    [x1,x2]=circlebc(t,1);   
    [dx1,dx2]=circlebc(t,2);
    [ddx1,ddx2]=circlebc(t,3);

else if bctype==2
    [x1,x2]=Penut(t,1);   
    [dx1,dx2]=Penut(t,2);
    [ddx1,ddx2]=Penut(t,3);
    else if bctype==3
            [x1,x2]=Pear(t,1);   
            [dx1,dx2]=Pear(t,2);
            [ddx1,ddx2]=Pear(t,3);
        else
            [x1,x2]=myrectangle(t,1);   
            [dx1,dx2]=myrectangle(t,2);
            [ddx1,ddx2]=myrectangle(t,3);
        end
    end
end
distance = sqrt( dx1.*dx1+dx2.*dx2 );

t_x1 = dx1./distance;
t_x2 = dx2./distance;
t_x = [t_x1', t_x2'];
n_x1 = t_x2;
n_x2 = -t_x1;
n_x =[n_x1', n_x2'];

%% scattering coefficient P_wave
alpha0=n_x*p0';
beta = kappa*alpha0 - sqrt(kappa^2*alpha0.^2-(kappa^2-1));

P11 = ones(n2,1)*p0(:,1)'-2*diag(n_x(:,1))*alpha0;
P12 = ones(n2,1)*p0(:,2)'-2*diag(n_x(:,2))*alpha0;

P21 = kappa*ones(n2,1)*p0(:,1)'- diag(n_x(:,1))*beta;
P22 = kappa*ones(n2,1)*p0(:,2)'- diag(n_x(:,2))*beta;

Pn21 = P22;
Pn22 = -P21;

Gamma1 = diag(n_x(:,1))*P11+diag(n_x(:,2))*P12;
Gamma2 = diag(n_x(:,1))*P21+diag(n_x(:,2))*P22;
Gamma3 = diag(n_x(:,1))*Pn21+diag(n_x(:,2))*Pn22;

A0 = 2*kappa*alpha0.^2-kappa-alpha0.*beta;
A1 = kappa - alpha0.*beta;
A2 = 2*alpha0.*(t_x*p0');

T01 = 1i*kp*(lambda*n_x(:,1)*ones(1,n1)+2*mu*alpha0*diag(p0(:,1)));
T02 = 1i*kp*(lambda*n_x(:,2)*ones(1,n1)+2*mu*alpha0*diag(p0(:,2)));

T11 =1i*kp*(lambda*n_x(:,1)*ones(1,n1)+2*mu*(Gamma1.*P11)).*A1./A0;
T12 =1i*kp*(lambda*n_x(:,2)*ones(1,n1)+2*mu*(Gamma1.*P12)).*A1./A0;

T21 = 1i*ks*mu*(Gamma2.*Pn21+Gamma3.*P21).*A2./A0;
T22 = 1i*ks*mu*(Gamma2.*Pn22+Gamma3.*P22).*A2./A0;

Rhp1 = real((T01+T11+T21)/(1i*kp));
Rhp2 = real((T02+T12+T22)/(1i*kp));

%% scattering coefficient S_wave
gamma = kappa1*alpha0 - sqrt(kappa1^2*alpha0.^2-(kappa1^2-1));
q0 = p0;
qn0 = pn0;

Q11 = kappa1*ones(n2,1)*q0(:,1)'-diag(n_x(:,1))*gamma;
Q12 = kappa1*ones(n2,1)*q0(:,2)'-diag(n_x(:,2))*gamma;
Q21 = ones(n2,1)*q0(:,1)'-2*diag(n_x(:,1))*alpha0;
Q22 = ones(n2,1)*q0(:,2)'-2*diag(n_x(:,2))*alpha0;
Qn21 = Q22;
Qn22 = -Q21;

Gamma1 = diag(n_x(:,1))*Q11+diag(n_x(:,2))*Q12;
Gamma2 = diag(n_x(:,1))*Q21+diag(n_x(:,2))*Q22;
Gamma3 = diag(n_x(:,1))*Qn21+diag(n_x(:,2))*Qn22;

B0 = 2*kappa1*alpha0.^2-kappa1-alpha0.*gamma;
B1 = -2*alpha0.*(t_x*q0');
B2 = kappa1 - alpha0.*gamma;
T01 =  1i*ks*mu*(alpha0*diag(qn0(:,1))+(n_x*qn0')*diag(q0(:,1)));
T02 =  1i*ks*mu*(alpha0*diag(qn0(:,2))+(n_x*qn0')*diag(q0(:,2))); 

T11 =  1i*kp*(lambda*n_x(:,1)*ones(1,n1)+2*mu*Gamma1.*Q11).*B1./B0;
T12 =  1i*kp*(lambda*n_x(:,2)*ones(1,n1)+2*mu*Gamma1.*Q12).*B1./B0;

T21 = 1i*ks*mu*(Gamma2.*Qn21+Gamma3.*Q21).*B2./B0;
T22 = 1i*ks*mu*(Gamma2.*Qn22+Gamma3.*Q22).*B2./B0;
    
Rhs1 = real((T01+T11+T21)/(1i*ks));
Rhs2 = real((T02+T12+T22)/(1i*ks));

for ni=1:n2
    for nj=1:n1
        if alpha0(ni,nj)>=0
            Rhp1(ni,nj)=0;
            Rhp2(ni,nj)=0;
            
            Rhs1(ni,nj)=0;
            Rhs2(ni,nj)=0;
        end
    end
end


 A = singlelayar_fullspace(omega, lambda, mu, n2/2, bctype);

 fp1 = exp(1i*kp*(x1'*cos(theta)+x2'*sin(theta)))*diag(cos(theta));
 fp2 = exp(1i*kp*(x1'*cos(theta)+x2'*sin(theta)))*diag(sin(theta));
 fp = [fp1;fp2];
 fs1 = exp(1i*ks*(x1'*cos(theta)+x2'*sin(theta)))*diag(sin(theta));
 fs2 =-exp(1i*ks*(x1'*cos(theta)+x2'*sin(theta)))*diag(cos(theta));
 fs = [fs1;fs2];

 
 f = [fp,fs];
 
 Phi = A\f;
 
 Rp1 = real(Phi(1:n2,1:n1)./(1i*kp*exp(1i*kp*(x1'*cos(theta)+x2'*sin(theta)))));
 Rp2 = real(Phi(n2+1:end,1:n1)./(1i*kp*exp(1i*kp*(x1'*cos(theta)+x2'*sin(theta)))));
 
 Rs1 = real(Phi(1:n2,n1+1:end)./(1i*ks*exp(1i*ks*(x1'*cos(theta)+x2'*sin(theta)))));
 Rs2 = real(Phi(n2+1:end,n1+1:end)./(1i*ks*exp(1i*ks*(x1'*cos(theta)+x2'*sin(theta)))));
%% p1
 figure, 
subplot(2,2,1);imagesc(t/pi,theta/pi,Rhp1');colorbar; 
colormap('jet');title(['Kirchhoff Approximation,omega=16pi']);
xlabel('Polar Coordinates of the Boundary, (Pi)'); ylabel('Incident Angle, (Pi)');

subplot(2,2,2);imagesc(t/pi,theta/pi,Rp1'); colorbar;
colormap('jet');title(['Real Scattering coefficienti,omega=16pi']);
xlabel('Polar Coordinates of the Boundary, (Pi)'); ylabel('Incident Angle, (Pi)');

subplot(2,2,3);
plot(t/pi,Rhp1(:,1),'-ob',t/pi,Rp1(:,1),'-*r'); 
legend('Kirchhoff','Real',4);
title('Incident Angle=0');
xlabel('Polar Coordinates of the Boundary, (Pi)');

subplot(2,2,4);
plot(t/pi,Rhp1(:,n1/2),'-ob',t/pi,Rp1(:,n1/2),'-*r'); 
legend('Kirchhoff','Real',4);
title('Incident Angle=Pi');
xlabel('Polar Coordinates of the Boundary, (Pi)');
 
 
%% p2
figure, 
subplot(2,2,1);imagesc(t/pi,theta/pi,Rhp2');colorbar; 
colormap('jet');title(['Kirchhoff Approximation,omega=16pi']);
xlabel('Polar Coordinates of the Boundary, (Pi)'); ylabel('Incident Angle, (Pi)');

subplot(2,2,2);imagesc(t/pi,theta/pi,Rp2'); colorbar;
colormap('jet');title(['Real Scattering coefficienti,omega=16pi']);
xlabel('Polar Coordinates of the Boundary, (Pi)'); ylabel('Incident Angle, (Pi)');

subplot(2,2,3);
plot(t/pi,Rhp2(:,1),'-ob',t/pi,Rp2(:,1),'-*r'); 
legend('Kirchhoff','Real',4);
title('Incident Angle=0');
xlabel('Polar Coordinates of the Boundary, (Pi)');

subplot(2,2,4);
plot(t/pi,Rhp2(:,n1/2),'-ob',t/pi,Rp2(:,n1/2),'-*r'); 
legend('Kirchhoff','Real',4);
title('Incident Angle=Pi');
xlabel('Polar Coordinates of the Boundary, (Pi)');
 
 %% s1
figure, 
subplot(2,2,1);imagesc(t/pi,theta/pi,Rhs1');colorbar; 
colormap('jet');title(['Kirchhoff Approximation,omega=16pi']);
xlabel('Polar Coordinates of the Boundary, (Pi)'); ylabel('Incident Angle, (Pi)');

subplot(2,2,2);imagesc(t/pi,theta/pi,Rs1'); colorbar;
colormap('jet');title(['Real Scattering coefficienti,omega=16pi']);
xlabel('Polar Coordinates of the Boundary, (Pi)'); ylabel('Incident Angle, (Pi)');

subplot(2,2,3);
plot(t/pi,Rhs1(:,1),'-ob',t/pi,Rs1(:,1),'-*r'); 
legend('Kirchhoff','Real',4);
title('Incident Angle=0');
xlabel('Polar Coordinates of the Boundary, (Pi)');

subplot(2,2,4);
plot(t/pi,Rhs1(:,n1/2),'-ob',t/pi,Rs1(:,n1/2),'-*r'); 
legend('Kirchhoff','Real',4);
title('Incident Angle=Pi');
xlabel('Polar Coordinates of the Boundary, (Pi)');

%% s2
figure, 
subplot(2,2,1);imagesc(t/pi,theta/pi,Rhs2');colorbar; 
colormap('jet');title(['Kirchhoff Approximation,omega=16pi']);
xlabel('Polar Coordinates of the Boundary, (Pi)'); ylabel('Incident Angle, (Pi)');

subplot(2,2,2);imagesc(t/pi,theta/pi,Rs2'); colorbar;
colormap('jet');title(['Real Scattering coefficienti,omega=16pi']);
xlabel('Polar Coordinates of the Boundary, (Pi)'); ylabel('Incident Angle, (Pi)');

subplot(2,2,3);
plot(t/pi,Rhs2(:,1),'-ob',t/pi,Rs2(:,1),'-*r'); 
legend('Kirchhoff','Real',4);
title('Incident Angle=0');
xlabel('Polar Coordinates of the Boundary, (Pi)');

subplot(2,2,4);
plot(t/pi,Rhs2(:,n1/2),'-ob',t/pi,Rs2(:,n1/2),'-*r'); 
legend('Kirchhoff','Real',4);
title('Incident Angle=Pi');
xlabel('Polar Coordinates of the Boundary, (Pi)');