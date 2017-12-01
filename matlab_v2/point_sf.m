%% ---------------- Parameters Setting---------------------------------------------%%


tic
n_recv =401;
lamda=1/2;
mu=1/4;
omega= 2*pi;
ks=omega/sqrt(mu);
kp=omega/sqrt(lamda+2*mu);

testpoint=[0;10];
receiver = zeros(2,n_recv);

a = 200;

receiver(1,:) = linspace(-a,a,n_recv);    receiver(2,:) =  zeros(1,n_recv);
   
Nx = 201;
Nz = 201;
NN = Nx*Nz;
x = linspace( -2, 2, Nx);
z = linspace( 8,12 , Nz);

G11=zeros(Nx,Nz);
G21=zeros(Nx,Nz);
G12=zeros(Nx,Nz);
G22=zeros(Nx,Nz);
%load pointdata G;
%% point scarttering data%%
G=Elastic_GreenTensor_Thalf_SIP5(omega,kp,ks,[testpoint(1)*ones(1,n_recv); testpoint(2)*ones(1,n_recv)],receiver);
G1=repmat(G(1,:),Nz,1);
G2=repmat(G(2,:),Nz,1);
G3=repmat(G(3,:),Nz,1);
G4=repmat(G(4,:),Nz,1);
%% 回传点源数据 %%
scaling = (receiver(1,end)-receiver(1,1))/n_recv;
%% 由于数据对称性，只需要计算一半
nx=(Nx-1)/2+1;
xh=x(1:nx);

g11=zeros(nx,Nz);
g21=zeros(nx,Nz);
g12=zeros(nx,Nz);
g22=zeros(nx,Nz);
 parfor jj=1:Nz  

           %gr=Traction_Dhalf(omega,kp,ks,[x(ix)*ones(1,n_recv); z(iz)*ones(1,n_recv)],receiver);
           re=[xh;repmat(z(jj),1,nx)];
           sample = reshape(repmat(re,n_recv,1),2,n_recv*nx);
           %Gr=TractionGreenTensor_2D(omega,kp,ks,sample,repmat(receiver,1,Nx));
           Gr=TractionDGreenTensor_2D_prin(omega,kp,ks,sample,repmat(receiver,1,nx));

           gr = zeros(nx,n_recv,4);
           for j=1:4 
               gr(:,:,j) = reshape(Gr(j,:),n_recv,nx).';
           end
           
           %gr=Elastic_GreenTensor_2D(omega,kp,ks,[x(ix)*ones(1,n_recv); z(iz)*ones(1,n_recv)],receiver);
           V1 = gr(:,:,1)*conj(G1(jj,:).')+gr(:,:,2)*conj(G2(jj,:).');
           V2 = gr(:,:,3)*conj(G1(jj,:).')+gr(:,:,4)*conj(G2(jj,:).');
           V3 = gr(:,:,1)*conj(G3(jj,:).')+gr(:,:,2)*conj(G4(jj,:).');
           V4 = gr(:,:,3)*conj(G3(jj,:).')+gr(:,:,4)*conj(G4(jj,:).');
           g11(:,jj)=scaling*V1;
           g21(:,jj)=scaling*V2;
           g12(:,jj)=scaling*V3;
           g22(:,jj)=scaling*V4;
 end

%% 利用对称性恢复矩阵
G11(1:nx,:)=g11;
G21(1:nx,:)=g21;
G12(1:nx,:)=g12;
G22(1:nx,:)=g22;

G11((nx+1):Nx,:)=G11((nx-1):-1:1,:);
G21((nx+1):Nx,:)=-G21((nx-1):-1:1,:);
G12((nx+1):Nx,:)=-G12((nx-1):-1:1,:);
G22((nx+1):Nx,:)=G22((nx-1):-1:1,:);

        
%% plot %%
[xx,zz]=meshgrid(x,z);
%save pointsf G11 G22 G12 G21

%figure, imagesc(x,z,real(G11)');colorbar;
figure
subplot(2,2,1)
 imagesc(x,z,-imag(G11)');colorbar;
colormap(jet)

subplot(2,2,2)
 imagesc(x,z,-imag(G12)');colorbar;
colormap(jet)

subplot(2,2,3)
 imagesc(x,z,-imag(G21)');colorbar;
colormap(jet)

subplot(2,2,4)
 imagesc(x,z,-imag(G22)');colorbar;
colormap(jet)


%figure, mesh(xx,zz,real(G11)');
%figure, mesh(xx,zz,-imag(G11)');
toc