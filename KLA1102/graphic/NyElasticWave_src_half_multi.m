function Data = NyElasticWave_src_half_multi(omegas, lambda, mu, n, bctype1,bctype2, n_src, n_recv, source, receiver)
%% the exact solution is given by the radiation solution
%% 2017 12 08 by zhou

w = quad_weights(n);
R = zeros(2*n);

for k=1:2*n
     idx=k:2*n;
     R(idx,k)=w(1:2*n-k+1);
     R(k,k)=R(k,k)/2;  %% for convinience
end
R=(R+R');

%% discrete point
node = 0:2*n-1;
t = pi*node(:)/n;
%图形1
if bctype1==1
    [x1,x2]=circlebcl(t,1);   
    [dx1,dx2]=circlebcl(t,2);

else if bctype1==2
    [x1,x2]=Penutl(t,1);   
    [dx1,dx2]=Penutl(t,2);
    else if bctype1==3
            [x1,x2]=p_leafl(t,1);   
            [dx1,dx2]=p_leafl(t,2);
        else
            [x1,x2]=myrectanglel(t,1);   
            [dx1,dx2]=myrectanglel(t,2);
        end
    end
end
%图形2
if bctype2==1
    [y1,y2]=circlebcr(t,1);   
    [dy1,dy2]=circlebcr(t,2);

else if bctype2==2
    [y1,y2]=Penutr(t,1);   
    [dy1,dy2]=Penutr(t,2);
    else if bctype2==3
            [y1,y2]=p_leafr(t,1);   
            [dy1,dy2]=p_leafr(t,2);
        else if bctype2==4
            [y1,y2]=myrectangler(t,1);   
            [dy1,dy2]=myrectangler(t,2);
            else
            [y1,y2]=kiter(t,1);   
            [dy1,dy2]=kiter(t,2);
            end
        end
    end
end


Ceuler = 0.577215664901532860;

distance = sqrt( dx1.*dx1+dx2.*dx2 );
distancey = sqrt( dy1.*dy1+dy2.*dy2 );

Nomegas = length(omegas);
Data = zeros(n_recv,n_src,2,2,Nomegas);
%% 对多频进行循环
for kk=1:Nomegas
    omega = omegas(kk);
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
    M{k}=zeros(4*n);
    H{k}=zeros(4*n);
    end
    tic
    %% 计算矩阵green tensor
    % 将矩阵拉开成向量
    %1,1块
    xl = repmat([x1,x2]',1,2*n);
    yl = reshape( repmat([x1,x2]',2*n,1),2,4*n*n);
    rGG  = Elastic_GreenTensor_Thalf_SIP_NY(omega,kp,ks,yl,xl);
    %2,2块
    xly = repmat([y1,y2]',1,2*n);
    yly = reshape( repmat([y1,y2]',2*n,1),2,4*n*n);
    rGGy  = Elastic_GreenTensor_Thalf_SIP_NY(omega,kp,ks,yly,xly);
    %1,2块
    rGG12  = Elastic_GreenTensor_Thalf_SIP_NY(omega,kp,ks,yly,xl);
    %2,1块
    rGG21  = Elastic_GreenTensor_Thalf_SIP_NY(omega,kp,ks,yl,xly);
    


    for j=1:2*n
        for k=1:2*n
            if (j==k)
                dist = distance(k);
                disty = distancey(k);
                
                g = rGG(:,(k-1)*2*n+j);
                G11 = reshape(g,2,2); 
                
                gy = rGGy(:,(k-1)*2*n+j);
                G22 = reshape(gy,2,2); 
                
                G1 = rGG12(:,(k-1)*2*n+j);
                G12 = reshape(G1,2,2);
                
                G2 = rGG21(:,(k-1)*2*n+j);
                G21 = reshape(G2,2,2);
            
                M{1,1}(j,k) = alpha; M{2,2}(j,k) = alpha;
                M{1,2}(j,k) = 0; M{2,1}(j,k) = 0;
                
                M{1,1}(j+2*n,k+2*n) = alpha; M{2,2}(j+2*n,k+2*n) = alpha;
                M{1,2}(j+2*n,k+2*n) = 0; M{2,1}(j+2*n,k+2*n) = 0;
            
                H{1,1}(j,k) = 2*alpha*log(dist) + ka1 + ka2*dx1(j)*dx1(j)/dist^2;
                H{2,2}(j,k) = 2*alpha*log(dist) + ka1 + ka2*dx2(j)*dx2(j)/dist^2;
            
                H{1,2}(j,k) =                           ka2*dx1(j)*dx2(j)/dist^2;
                H{2,1}(j,k) =                           ka2*dx2(j)*dx1(j)/dist^2;
                
                H{1,1}(j+2*n,k+2*n) = 2*alpha*log(disty) + ka1 + ka2*dy1(j)*dy1(j)/disty^2;
                H{2,2}(j+2*n,k+2*n) = 2*alpha*log(disty) + ka1 + ka2*dy2(j)*dy2(j)/disty^2;
            
                H{1,2}(j+2*n,k+2*n) =                           ka2*dy1(j)*dy2(j)/disty^2;
                H{2,1}(j+2*n,k+2*n) =                           ka2*dy2(j)*dy1(j)/disty^2;
                for l=1:2
                    for m=1:2
                        M{l,m}(j,k) = dist*M{l,m}(j,k);
                        H{l,m}(j,k) = dist*(H{l,m}(j,k)++G11(l,m));
                        
                        H{l,m}(j,k+2*n) = disty*G12(l,m);
                        
                        M{l,m}(j+2*n,k+2*n) = disty*M{l,m}(j+2*n,k+2*n);
                        H{l,m}(j+2*n,k+2*n) = disty*(H{l,m}(j+2*n,k+2*n)+G22(l,m));
                        
                        H{l,m}(j+2*n,k) = dist*G21(l,m);
                    end
                end
            else
                dist = distance(k);
                disty = distancey(k);
                v = ([x1(j)-x1(k) x2(j)-x2(k)]); 
                vy = ([y1(j)-y1(k) y2(j)-y2(k)]);
                lg4s = log(4*sin((t(j)-t(k))/2)^2) ; 
                rv = norm(v);
                rvy = norm(vy);
                
                g = rGG(:,(k-1)*2*n+j);
                G11 = reshape(g,2,2); 
                
                gy = rGGy(:,(k-1)*2*n+j);
                G22 = reshape(gy,2,2); 
                
                G1 = rGG12(:,(k-1)*2*n+j);
                G12 = reshape(G1,2,2);
                
                G2 = rGG21(:,(k-1)*2*n+j);
                G21 = reshape(G2,2,2);
            
                phi1 = 1/(2*pi)*( -1/(2*mu)*besselj(0,ks*rv) + 1/(2*w2*rv)* ( ks*besselj(1,ks*rv) - kp*besselj(1,kp*rv) ) ) ;
                phi2 = 1/(2*pi)*( 1/(2*w2)*( ks^2*besselj(0,ks*rv) -2*ks/rv*besselj(1,ks*rv) - kp^2*besselj(0,kp*rv) + 2*kp/rv*besselj(1,kp*rv)  ) );
                
                phiy1 = 1/(2*pi)*( -1/(2*mu)*besselj(0,ks*rvy) + 1/(2*w2*rvy)* ( ks*besselj(1,ks*rvy) - kp*besselj(1,kp*rvy) ) ) ;
                phiy2 = 1/(2*pi)*( 1/(2*w2)*( ks^2*besselj(0,ks*rvy) -2*ks/rvy*besselj(1,ks*rvy) - kp^2*besselj(0,kp*rvy) + 2*kp/rvy*besselj(1,kp*rvy)  ) );
                for l=1:2
                    for m=1:2
                        e=0;
                        if l==m
                            e=1;
                        end
                    
                        M{l,m}(j,k) = dist*(phi1*e + phi2*v(l)*v(m)/rv^2);
                        H{l,m}(j,k) = (G11(l,m)*dist -   M{l,m}(j,k)*lg4s);
                        
                        M{l,m}(j+2*n,k+2*n) = disty*(phiy1*e + phiy2*vy(l)*vy(m)/rvy^2);
                        H{l,m}(j+2*n,k+2*n) = (G22(l,m)*disty -   M{l,m}(j+2*n,k+2*n)*lg4s);
                        
                        H{l,m}(j,k+2*n) = disty*G12(l,m);
                        H{l,m}(j+2*n,k) = dist*G21(l,m);
                    end
                end
           
            end
        end
    end

    toc
    disp('参数矩阵生成');
    %% 构造线性系统，方程左端的矩阵 %%
    
    A = zeros(8*n,8*n);

    A(1:4*n,1:4*n)  = repmat(R,2,2).*M{1,1} + pi/n*H{1,1};
    A(1:4*n,4*n+1:end) = repmat(R,2,2).*M{1,2} + pi/n*H{1,2};
    A(4*n+1:end,1:4*n) = repmat(R,2,2).*M{2,1} + pi/n*H{2,1};
    A(4*n+1:end,4*n+1:end) = repmat(R,2,2).*M{2,2} + pi/n*H{2,2};
    tic
    %% 计算线性方程右端项
    f=zeros(8*n,2*n_src);
    % 同样的将循环去掉，用长向量代替
    xsource = reshape(repmat(source(:,:),2*n,1),2,2*n*n_src);
    xobs = repmat([x1,x2]',1,n_src);
    yobs = repmat([y1,y2]',1,n_src);
    % x,y 互换 转置
    Gs1 = Elastic_GreenTensor_Thalf_SIP4(omega,kp,ks,xobs,xsource);
    Gsy1 = Elastic_GreenTensor_Thalf_SIP4(omega,kp,ks,yobs,xsource);
    Gs = - transpose(Gs1); 
    Gsy = - transpose(Gsy1); 
    % 右端项为点源两个方向的组合，每一列代表一个方向
    f(1:2*n,:) = [reshape(Gs(:,1),2*n ,n_src),reshape(Gs(:,2),2*n ,n_src)];
    f(2*n+1:4*n,:)=[reshape(Gsy(:,1),2*n ,n_src),reshape(Gsy(:,2),2*n ,n_src)];
    f(4*n+1:6*n,:) = [reshape(Gs(:,3), 2*n ,n_src),reshape(Gs(:,4), 2*n ,n_src)];
    f(6*n+1:end,:) = [reshape(Gsy(:,3), 2*n ,n_src),reshape(Gsy(:,4), 2*n ,n_src)];
    toc
    disp('右端项计算完成');
    
    
    
    %% 解由单层位势表示的线性系统  
    %% U = S\phi
    
    
    % Phi是包含两个方向点源产生的核，前n_src代表一个方向，其余代表一个方向
    tic
    Phi = A\f;
    
    toc
    disp('得到线性方程解，即位势核');
    
   %% 计算样本点的散射数据
    % 同样的将循环去掉，用长向量代替
    tic
    if n_src==n_recv
        
        Gr = Gs1;
        Gry = Gsy1;

    else
        xreceiver = reshape(repmat(receiver(:,:),2*n,1),2,2*n*n_recv);
        xobs2= repmat([x1,x2]',1,n_recv);
        yobs2= repmat([y1,y2]',1,n_recv);
        Gr = Elastic_GreenTensor_Thalf_SIP4(omega,kp,ks,xobs2,xreceiver);
        Gry = Elastic_GreenTensor_Thalf_SIP4(omega,kp,ks,yobs2,xreceiver);
    end
    
    G1=reshape(Gr(1,:),2*n,n_recv).';
    G2=reshape(Gr(2,:),2*n,n_recv).';
    G3=reshape(Gr(3,:),2*n,n_recv).';
    G4=reshape(Gr(4,:),2*n,n_recv).';
    distances = repmat(distance(:),1,2*n_src);
    
    Gy1=reshape(Gry(1,:),2*n,n_recv).';
    Gy2=reshape(Gry(2,:),2*n,n_recv).';
    Gy3=reshape(Gry(3,:),2*n,n_recv).';
    Gy4=reshape(Gry(4,:),2*n,n_recv).';
    distancesy = repmat(distancey(:),1,2*n_src);
    
    toc
    disp('单层位势函数（障碍物到接收点）计算完成');
    
    tic
    
    temp1 = G1*(distances.*Phi(1:2*n,:)) + G3*(distances.*Phi(4*n+1:6*n,:));
    temp2 = G2*(distances.*Phi(1:2*n,:)) + G4*(distances.*Phi(4*n+1:6*n,:));
    
    tempy1 = Gy1*(distancesy.*Phi(2*n+1:4*n,:)) + Gy3*(distances.*Phi(6*n+1:end,:));
    tempy2 = Gy2*(distancesy.*Phi(2*n+1:4*n,:)) + Gy4*(distances.*Phi(6*n+1:end,:));
    
    
    
        Data(:,:,1,1,kk)=pi/n*(temp1(:,1:n_src)+tempy1(:,1:n_src));
        Data(:,:,2,1,kk)=pi/n*(temp2(:,1:n_src)+tempy2(:,1:n_src));
        Data(:,:,1,2,kk)=pi/n*(temp1(:,(n_src+1):2*n_src)+tempy1(:,(n_src+1):2*n_src));
        Data(:,:,2,2,kk)=pi/n*(temp2(:,(n_src+1):2*n_src)+tempy2(:,(n_src+1):2*n_src));
    
    toc
end

disp('数据合成结束');
end