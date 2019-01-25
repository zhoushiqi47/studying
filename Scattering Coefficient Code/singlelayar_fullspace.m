function A=singlelayar_fullspace(omega, lambda, mu, n, bctype)

w = quad_weights(n);
R = zeros(2*n);
lamda=lambda;
for k=1:2*n
     idx=k:2*n;
     R(idx,k)=w(1:2*n-k+1);
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
    [x1,x2]=Penut(t,1);   
    [dx1,dx2]=Penut(t,2);
    else if bctype==3
            [x1,x2]=Pear(t,1);   
            [dx1,dx2]=Pear(t,2);
        else
            [x1,x2]=myrectangle(t,1);   
            [dx1,dx2]=myrectangle(t,2);
        end
    end
end

Ceuler = 0.577215665;
distance = sqrt( dx1.*dx1+dx2.*dx2 );

    cp = sqrt( lambda + 2*mu);
    cs = sqrt( mu );
    %% the corresponding wavenumber;
    kp = omega / cp;
    ks = omega / cs;

    %% constant definition
    w2 = omega^2;
    alpha = 1/(2*pi)*( - 1/(4*w2)*(ks^2 + kp^2) );
    % beta1 = 1/(2*pi)*( 1/(32*w2)*( 3*ks^4 + kp^4 ) );
    % beta2 = 1/(2*pi)*( 1/(16*w2)*(   kp^4 - ks^4 ) );
    ka1 = -1/(4*pi*w2)*( ks^2*log(ks/2) + kp^2*log(kp/2) + 1/2*(ks^2 - kp^2) + (Ceuler - 1i*pi/2)*(ks^2+kp^2)  );
    ka2 = 1/(4*pi*w2)*( ks^2 - kp^2 );
    %% end of constant definiation

    M = cell(2,2);
    H = cell(2,2);
    for k=1:4
    M{k}=zeros(2*n);
    H{k}=zeros(2*n);
    end

    %% ????????green tensor
    % ????????????????
    xl = repmat([x1,x2]',1,2*n);
    yl = reshape( repmat([x1,x2]',2*n,1),2,4*n*n);
    %rGG  = Elastic_GreenTensor_Thalf_SIP_NY(omega,kp,ks,yl,xl);
    rGG = Elastic_GreenTensor_2D(omega,kp,ks,xl,yl);


    
    
    
    
    for j=1:2*n
        for k=1:2*n
            if (j==k)
                dist = distance(k);

                M{1,1}(j,k) = alpha; M{2,2}(j,k) = alpha;
                M{1,2}(j,k) = 0; M{2,1}(j,k) = 0;
            
                H{1,1}(j,k) = 2*alpha*log(dist) + ka1 + ka2*dx1(j)*dx1(j)/dist^2;
                H{2,2}(j,k) = 2*alpha*log(dist) + ka1 + ka2*dx2(j)*dx2(j)/dist^2;
            
                H{1,2}(j,k) =                           ka2*dx1(j)*dx2(j)/dist^2;
                H{2,1}(j,k) =                           ka2*dx2(j)*dx1(j)/dist^2;
                for l=1:2
                    for m=1:2
                        M{l,m}(j,k) = dist*M{l,m}(j,k);
                        H{l,m}(j,k) = dist*H{l,m}(j,k);
                    end
                end
            else
                dist = distance(k);
                v = ([x1(j)-x1(k) x2(j)-x2(k)]); 
                lg4s = log(4*sin((t(j)-t(k))/2)^2) ; 
                rv = norm(v);
                G1 = rGG(:,(k-1)*2*n+j);
                G = reshape(G1,2,2); 
            
                phi1 = 1/(2*pi)*( -1/(2*mu)*besselj(0,ks*rv) + 1/(2*w2*rv)* ( ks*besselj(1,ks*rv) - kp*besselj(1,kp*rv) ) ) ;
                phi2 = 1/(2*pi)*( 1/(2*w2)*( ks^2*besselj(0,ks*rv) -2*ks/rv*besselj(1,ks*rv) - kp^2*besselj(0,kp*rv) + 2*kp/rv*besselj(1,kp*rv)  ) );
                for l=1:2
                    for m=1:2
                        e=0;
                        if l==m
                            e=1;
                        end
                    
                        M{l,m}(j,k) = dist*(phi1*e + phi2*v(l)*v(m)/rv^2);
                        H{l,m}(j,k) = (G(l,m)*dist -   M{l,m}(j,k)*lg4s);
                    end
                end
           
            end
        end
    end



    %% ???????????????????????????? %%
    
    A = zeros(4*n,4*n);

    A(1:2*n,1:2*n)  = R.*M{1,1} + pi/n*H{1,1};
    A(1:2*n,2*n+1:end) = R.*M{1,2} + pi/n*H{1,2};
    A(2*n+1:end,1:2*n) = R.*M{2,1} + pi/n*H{2,1};
    A(2*n+1:end,2*n+1:end) = R.*M{2,2} + pi/n*H{2,2};
end
