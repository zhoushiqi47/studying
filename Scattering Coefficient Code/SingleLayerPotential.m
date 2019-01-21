function [Mat,x1,x2,dx1,dx2,distance]=SingleLayerPotential(n, bctype,wavenumber)
%%  source: positon of source
%% receive: position of receive

%% compute the quadrature weights;

w = quad_weights(n);
R = zeros(2*n);

for k=1:2*n
     idx=[k:2*n];
     R(idx,k)=w([1:2*n-k+1]);
     R(k,k)=R(k,k)/2;  %% for convinience
end
R=(R+R');
%% discrete point
node = 0:2*n-1;
t = pi*node(:)/n;
if bctype==1
    [x1,x2]=circlebc(t,1);   
    [dx1,dx2]=circlebc(t,2);

else if bctype==2
    [x1,x2]=kite(t,1);   
    [dx1,dx2]=kite(t,2);
    else if bctype==3
               [x1,x2]=leaf(t,1);   
               [dx1,dx2]=leaf(t,2);
        else if bctype==4
            [x1,x2]=Pear(t,1);   
            [dx1,dx2]=Pear(t,2);
            else
                [x1,x2]=rectanglebc(t,1);   
                [dx1,dx2]=rectanglebc(t,2);
            end
        end
    end

end

r = zeros(2*n);

for k=1:2*n
    for j=1:2*n
        r(k,j)=sqrt((x1(j)-x1(k))* (x1(j)-x1(k))+(x2(j)-x2(k))* (x2(j)-x2(k)));
    end
end


M1 = zeros(2*n);
M2 = zeros(2*n);
Ceuler = 0.577215665;
distance = sqrt( dx1.*dx1+dx2.*dx2 );

for j=1:2*n
    for k=1:2*n
        if (j==k)
            dist = distance(k);
            M1(j,k) = -1/(2*pi)*besselj(0,wavenumber*r(j,k))*dist;
            M2(j,k) = (i/2-Ceuler/pi-1/(2*pi)*log(wavenumber^2/4*dist*dist))*dist;
        else
            dist = distance(k);
            temp = dx2(k)*(x1(j)-x1(k))-dx1(k)*(x2(j)-x2(k));
            M1(j,k) = -1/(2*pi)*besselj(0,wavenumber*r(j,k))*dist;
            M2(j,k) = i/2*besselh(0,wavenumber*r(j,k))*dist - M1(j,k)*log(4*sin((t(j)-t(k))/2)^2);
        end
    end
end
%% the linear System is A = I-(L1+ik*M1 +pi/n*(L2+ik*M2))

%% Mat = 2*(R.*M1 + pi/n*M2);
Mat = R.*M1 + pi/n*M2;
% eta = wavenumber; L1=0;L2=0;
% A=eye(2*n) - (R.*(L1+1i*eta*M1)+pi/n*(L2+1i*eta*M2));
% A=-A;
%% the right hand side is double of incidient 

%% 
% for is=1:n_src
%     f(:,is) = -2*Green(wavenumber, source(:,is)*ones(1,2*n),[x1';x2']);
% end

%% the solution is the potential on the boundary of D
% phi = A\f;
% 
% %% Composite Trapzitol Formula for Computing Far fields pattern
% 
% 
% 
% 
% U = zeros(n_recv,n_src);
% for j=1:n_src
%     for k=1:n_recv
%        
%         g1 = Green(wavenumber,receiver(:,k)*ones(1,2*n),[x1';x2']);
%         gd1 = Green_Grad(wavenumber,receiver(:,k)*ones(1,2*n),[x1';x2']);
%         temp = ( (gd1(1,:).*dx2' - gd1(2,:).*dx1') - (1i* eta*distance').*g1)*phi(:,j);
%         U(k,j)=2*pi/n*temp;
%     end
% end
% Nx = 201;
% Nz = 201;
% Irtm = zeros(Nx,Nz);
% x = linspace( -8, 8, Nx);
% z = linspace( -8, 8, Nz);
% Nfreq = 1;size(U)
% for ix = 1 : Nx
%     for iz = 1: Nz
%         for k=1:Nfreq
%             gr = conj(((Green(wavenumber,receiver,[x(ix)*ones(1,n_recv); z(iz)*ones(1,n_recv)]))));
%             gs = conj(((Green(wavenumber,source,[x(ix)*ones(1,n_src); z(iz)*ones(1,n_src)]))));
%             Irtm(ix,iz)=Irtm(ix,iz)+sum( gs.* ( gr*U(:,:) ) );
%         end
%     end
% end
% 
% 
% figure
% [xx,zz]=meshgrid(x,z);
% subplot(1,3,1)
% imagesc(abs(imag(Irtm))');
% subplot(1,3,2);
% mesh(xx,zz,abs(imag(Irtm))');
% subplot(1,3,3)
% plot(x,imag(Irtm(101,:)));