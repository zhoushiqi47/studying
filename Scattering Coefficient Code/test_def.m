clear;
f=1;
freq=pi*f;
radius=1;
flag=1;
n1=128; %% discretization of incidence angle
gamma = 2*pi/n1*(0:n1-1);
eta = [cos(gamma); sin(gamma)]';
bctype = 1;

n2 = 256; %% discretization of the boundary of obstacle
theta = 2*pi/n2*(0:n2-1);
if bctype==1
    [z1,z2]=circlebc(theta,2);   
 
else if bctype ==2
    [z1,z2]=kite(theta,2);   
    else if bctype ==3
        [z1,z2]=leaf(theta,2);
        else if bctype==4
            [z1,z2]=Pear(theta,2);
            else
                [z1,z2]=rectanglebc(theta,2);
            end
        end
    end
end
z = sqrt(z1.^2+z2.^2);

y=[(z2./z).' (-z1./z).'];

%% Kirchhoff approximation of scattering coefficient
rt = 2*y*eta';
for ni=1:n2
    for nj=1:n1
        if rt(ni,nj)>0
            rt(ni,nj)=0;
        end
    end
end

%% By definition of scattering coefficient of sound soft non-penetrable obstacle 
npts = n2/2;
tic

[Mat,x1,x2,dx1,dx2,distance]=SingleLayerPotential(npts, bctype,freq);

    
[l,u]=lu(Mat);

F0 = exp(1i*freq*([x1(:) x2(:)]*eta'));
%% dudv = -u\(l\F0);
dudv = u\(l\(2*F0));

rtt = (-1i/freq)*dudv./F0;

rtt = real(rtt);
c1 = min(-2,min(min(rtt)));
c2 = max(0,max(max(rtt)));


figure,
subplot(2,2,1);imagesc(theta/pi,gamma/pi,rt');colorbar; 
caxis([c1,c2]);
colormap('jet');title(['Kirchhoff Approximation,k=', num2str(f),'pi']);
xlabel('Polar Coordinates of the Boundary, (Pi)'); ylabel('Incident Angle, (Pi)')

subplot(2,2,2);imagesc(theta/pi,gamma/pi,real(rtt)'); colorbar;
 caxis([c1,c2]);
colormap('jet');title(['Real Scattering coefficient,k=', num2str(f),'pi']);
xlabel('Polar Coordinates of the Boundary, (Pi)'); ylabel('Incident Angle, (Pi)')
subplot(2,2,3);
plot(theta/pi,rt(:,1),'-ob',theta/pi,rtt(:,1),'-*r'); 
legend('Kirchhoff','Real',4);
title('Incident Angle=0');
xlabel('Polar Coordinates of the Boundary, (Pi)');
subplot(2,2,4);
plot(theta/pi,rt(:,n1/2),'-ob',theta/pi,rtt(:,n1/2),'-*r'); 
legend('Kirchhoff','Real',4);
title('Incident Angle=Pi');
xlabel('Polar Coordinates of the Boundary, (Pi)');

%saveas(gcf,['scattering coefficient leaf_',num2str(f),'_pi.eps']);
% % Plot in other ways
% figure,imagesc(real(rtt)');colorbar;colormap('jet');
% figure,imagesc(imag(rtt)');colorbar;colormap('jet');
% 
% figure,
% subplot(1,2,1);imagesc(theta/pi,gamma/pi,rt');colorbar;colormap('jet');title('Kirchhoff Approximation');
% subplot(1,2,2);imagesc(theta/pi,gamma/pi,real(rtt)'); colorbar;colormap('jet');title('Computation by Us');
% 
% figure,
% subplot(1,2,1);
% plot(theta/pi,rt(:,1),'-ob');
% hold on
% plot(theta/pi,rtt(:,1),'-*r');
% hold off
% title('theta=0');
% subplot(1,2,2);
% plot(theta/pi,rt(:,n1/2),'-ob');
% hold on
% plot(theta/pi,rtt(:,n1/2),'-*r');
% hold off
% title('theta=pi');