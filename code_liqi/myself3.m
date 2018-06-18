function power=myself3(n)
h=0.21;
lambd=[0.1,0.1];
k=1000;
p1=1;
p2=2;
ss=0;
re=400;
for l=1:k
    fprintf('l=%d\n',l);
    alphabet = [0 1]; prob = [0.5 0.5];
    w=randsrc(n,p2,[alphabet; prob]);
    x=randn(n,p1);
    error=randn(n,1);  
    z=[x w];
%---------------假设--------------------
 beta0=[1 1 1 1 1]; 
 y=beta0(1)+beta0(2)*w(:,1)+beta0(3)*w(:,2)+beta0(4)*x(:,1)+beta0(5)*x(:,1).^2+error;    % DGP0
%  y=beta0(1)+beta0(2)*w(1,:)+beta0(3)*w(2,:)+beta0(4)*x(1,:)+beta0(5)*x(1,:).^2+0.25*x(1,:).^3+error; % DGP1
%  y=beta0(1)+beta0(2)*x(1,:)+beta0(3)*x(2,:)+beta0(4)*z(1,:)+beta0(5)*z(1,:).^2+0.5*x(1,:).*z(1,:)+0.5*x(2,:).*z(1,:).^2+error; % DGP2
%  y=beta0(1)+beta0(2)*x(1,:)+beta0(3)*x(2,:)+beta0(4)*z(1,:)+beta0(5)*z(1,:).^2+(x(1,:)+x(2,:)).*sin(4*pi*z(1,:))+error;  % DGP3

%----------------------------------------

%    估计beta
beta0=[1 1 1 1 1]./2; 
 beta=nlinfit(z,y,@myselfmodel,beta0);
 resi=y-(beta(1)+beta(2)*z(:,2)+beta(3)*z(:,3)+beta(4)*z(:,1)+beta(5)*z(:,1).^2);
%  [a,b]=CV(w,y,x,n,p2,p1);
 
 %      L函数
%   lambd=b*ones(p2,1);
% lambd=0*ones(n,1);%./10;
Lx=zeros(n,n*p2);
Lz=zeros(n,n*p2);
L=zeros(n,n*p2);
for s=1:p2
    Lz(:,(s-1)*n+1:s*n)=lambd(s).^abs(w(:,s)*ones(1,n)-ones(n,1)*w(:,s)');
end
L=Lz(:,1:n);
for s=2:p2
    L=L.*Lz(:,(s-1)*n+1:s*n);
end

%   W函数
% h=a*ones(p1,1);
% h=ones(m,1).*p^(-1/(4+m));
W=zeros(n,n*p1);
Wz=zeros(n,n*p1);
for s=1:p1
    Wx=(x(:,s)*ones(1,n)-ones(n,1)*x(:,s)')/h(s);
    Wz(:,(s-1)*n+1:s*n)=exp(-Wx.^2/2)/h(s);
end
W=Wz(:,1:n);
for s=2:p1
    W=W.*Wz(:,(s-1)*n+1:s*n);
end

%      L.*W
   KZ=L.*W;
   AZ=KZ.*(ones(n,1)*resi');
   AZ1=sum(AZ,2)-diag(AZ);
   Vn=resi'*AZ1/(n^2);
    
   BZ=KZ.*KZ.*(ones(n,1)*resi.^2');
   BZ1=sum(BZ,2)-diag(BZ);
   Sn=resi.^2'*BZ1;
   Sn=2*prod(h)*Sn/(n^2);
 
   TfinalZ=n*sqrt(prod(h))*Vn/sqrt(Sn);
   
%-----------bootstrap-----------
   for t=1:re
%        fprintf('t=%d\n',t);
       prob=(1+sqrt(5))/(2*sqrt(5));
       Z=random('binomial',1,(1-prob)*ones(n,1));
       a1=0.5*(1-sqrt(5));a2=0.5*(1+sqrt(5));
       V=(a2-a1)*Z+a1;
       err_star=resi.*V;
       y_star=beta(1)+beta(2)*z(:,2)+beta(3)*z(:,3)+beta(4)*z(:,1)+beta(5)*z(:,1).^2+err_star;
       beta_star=nlinfit(z,y_star,@myselfmodel,beta0);
% %        if (isnan(beta)==1)|(isinf(beta)==1)
% %            m=m+1;
% %            continue; 
% %        end
       resi_star=y_star-(beta_star(1)+beta_star(2)*z(:,2)+beta_star(3)*z(:,3)+beta_star(4)*z(:,1)+beta_star(5)*z(:,1).^2);

       AZ=KZ.*(ones(n,1)*resi_star');
       AZ1=sum(AZ,2)-diag(AZ);
       Vn=resi_star'*AZ1/(n^2);
       
       BZ=KZ.*KZ.*(ones(n,1)*resi_star.^2');
       BZ1=sum(BZ,2)-diag(BZ);
       Sn=resi_star.^2'*BZ1;
       Sn=2*prod(h)*Sn/(n^2);
       
       TstarZZ(t)=n*sqrt(prod(h))*Vn/sqrt(Sn);
   
   end     
   TstarZ=sort(TstarZZ);
%    critic=TstarZ(1,floor((19/20)*re));
   critic=TstarZ(1,(19/20)*re);
   powersZ(l)=(TfinalZ>critic);
   mean(powersZ)
end
power=mean(powersZ)
end