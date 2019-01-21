clear;
f=1;
freq=pi*f;
radius=1;
flag=1;
n1=128; %% discretization of incidence angle
gamma = 2*pi/n1*(0:n1-1);
eta = [cos(gamma); sin(gamma)]';
bctype = 3;

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
if flag==1
    [Mat,x1,x2,dx1,dx2,distance]=SingleLayerPotential(npts, bctype,freq);
else
    [Mat,x1,x2,dx1,dx2,distance]=CombineLayerpotential(npts, bctype,freq);
end
    
[l,u]=lu(Mat);

F0 = exp(1i*freq*([x1(:) x2(:)]*eta'));
%% dudv = -u\(l\F0);
dudv = u\(l\(2*F0));
if flag==1
    rtt = (-1i/freq)*dudv./F0;

else
    dus=zeros(npts*2,n1);
    dui=zeros(npts*2,n1);
    X=[x1(:)';x2(:)'];DX=[dx1(:)';dx2(:)'];
    Y=X;DY=DX;
    for id=1:n1
        gd2=zeros(npts*2,npts*2);
        gd1=zeros(npts*2,npts*2);
                d=sqrt(X'*Y);
                vad0= besselh(0,freq*d);
                vad1= besselh(1,freq*d);
                vad2= besselh(2,freq*d);
                temp=0.5*(vad0-vad2);
                if d==0
                    gd1=0;
                    gd2=0;                   
                else
                    gd1=-1i/4*vad1*X'*DX./d;
                    gd2=1i/4*temp.*freq^2*X'*DX*Y'*DY./d./d-1i/8*temp.*freq*X'*DX*Y'*DY./d./d./d./d;
                end    
             dus(:,id)=(gd2-1i*freq*gd1)*dudv(:,id)+1i/2*freq*dudv(:,id);
             dui(:,id)=exp(1i*freq*eta(id,:)*X)*1i*freq.*(eta(id,:)*DX);  
%         for ixhat=1:npts*2
%                 d=sqrt(X(:,ixhat)'*Y);
%                 vad0= besselh(0,freq*d);
%                 vad1= besselh(1,freq*d);
%                 vad2= besselh(2,freq*d);
%                 temp=0.5*(vad0-vad2);
%                 if d==0
%                     gd1(ixhat,:)=0;
%                     gd2(ixhat,:)=0;                   
%                 else
%                     gd1(ixhat,:)=-1i/4*vad1*X(:,ixhat)'*DX(:,ixhat)./d;
%                     gd2(ixhat,:)=1i/4*temp.*freq^2*X(:,ixhat)'*DX(:,ixhat)*Y'*DY./d./d-1i/8*temp.*freq*X(:,ixhat)'*DX(:,ixhat)*Y'*DY./d./d./d./d;
%                 end    
% %             for iy=1:npts*2
% %                 d=sqrt(X(:,ixhat)'*Y(:,iy));
% %                 vad0= besselh(0,freq*d);
% %                 vad1= besselh(1,freq*d);
% %                 vad2= besselh(2,freq*d);
% %                 temp=0.5*(vad0-vad2);
% %                 if d==0
% %                     gd1(ixhat,iy)=0;
% %                     gd2(ixhat,iy)=0;                   
% %                 else
% %                     gd1(ixhat,iy)=-1i/4*vad1*X(:,ixhat)'*DX(:,ixhat)/d;
% %                     gd2(ixhat,iy)=1i/4*temp*freq^2*X(:,ixhat)'*DX(:,ixhat)*Y(:,iy)'*DY(:,iy)/d/d-1i/8*temp*freq*X(:,ixhat)'*DX(:,ixhat)*Y(:,iy)'*DY(:,iy)/d/d/d/d;
% %                 end                
% %             end
%             dus(ixhat,id)=(gd2(ixhat,:)-1i*freq*gd1(ixhat,:))*dudv(:,id)+1i/2*freq*dudv(ixhat,id);
%             dui(ixhat,id)=exp(1i*freq*eta(id,:)*(X(:,ixhat)))*1i*freq*eta(id,:)*DX(:,ixhat);          
%         end
    end
    dudv=dus+dui;
    rtt = (-1i/freq)*dudv./F0;
   

end
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