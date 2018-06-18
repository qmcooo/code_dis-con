function myselfnull()
tic
k=1000;
n=200;
p=2;
q=1;
ss=0;
re=399;
for l=1:k
    fprintf('l=%d\n',l);
alphabet = [0 1]; prob = [0.5 0.5];
x=randsrc(p,n,[alphabet; prob]);
%x=randn(n,p);
% z=normrnd(0,1,m,p);
z=randn(q,n);
beta0=[1 1 1 1 1];
%  beta0=[1 1 1 1 1 1];
%     y    
error=randn(1,n);  
 y=beta0(1)+beta0(2)*x(1,:)+beta0(3)*x(2,:)+beta0(4)*z(1,:)+beta0(5)*z(1,:).^2+error;    % DGP0
% y=beta0(1)+beta0(2)*x(1,:)+beta0(3)*x(2,:)+beta0(4)*z(1,:)+beta0(5)*z(1,:).^2+beta0(6)*z(1,:).^3+error; % DGP1

%    ¹À¼Æbeta
xz=[x;z];
%fun=@(beta,xz)beta(1)+beta(2)*xz(1,:)+beta(3)*xz(2,:)+beta(4)*xz(3,:)+beta(5)*xz(3,:).^2;  
% fun=@(beta,xz)beta(1)+beta(2)*xz(1,:)+beta(3)*xz(2,:)+beta(4)*xz(3,:)+beta(5)*xz(3,:).^2+beta(6)*xz(3,:).^3;

beta=nlinfit(xz,y,@myselfmodel,beta0);
% beta=lsqcurvefit(fun, beta0, xz, y);

%     residual
 resi=y-(beta(1)+beta(2)*x(1,:)+beta(3)*x(2,:)+beta(4)*z+beta(5)*z.^2);
%  resi=y-(beta(1)+beta(2)*x(1,:)+beta(3)*x(2,:)+beta(4)*z+beta(5)*z.^2+beta(6)*z(1,:).^3);
 [a,b]=CV(x,y,z,n,p,q);
 
 %      Lº¯Êý
  lambd=b*ones(p,1);
% lambd=0*ones(n,1);%./10;
Lx=zeros(n,n*p);
Lz=zeros(n,n*p);
L=zeros(n,n*p);
for s=1:p
    Lx(:,(s-1)*n+1:s*n)=x(s,:)'*ones(1,n)-ones(n,1)*x(s,:);
    Lz(:,(s-1)*n+1:s*n)=lambd(s).^(1-(Lx(:,(s-1)*n+1:s*n)==0));
end
L=Lz(:,1:n);
for s=2:p
    L=L.*Lz(:,(s-1)*n+1:s*n);
end

%   Wº¯Êý
h=a*ones(q,1);
% h=ones(m,1).*p^(-1/(4+m));
W=zeros(n,n*q);
Wz=zeros(n,n*q);
for s=1:q
    Wx(:,(s-1)*n+1:s*n)=z(s,:)'*ones(1,n)-ones(n,1)*z(s,:);
    Wz(:,(s-1)*n+1:s*n)=exp(-((Wx(:,(s-1)*n+1:s*n)/h(s)).^2)./2)/h(s);
end
W=Wz(:,1:n);
for s=2:q
    W=W.*Wz(:,(s-1)*n+1:s*n);
end

%      L.*W
   KZ=L.*W;
   AZ=KZ.*(ones(n,1)*resi);
   AZ1=sum(AZ,2)-diag(AZ);
   AZ2=resi*AZ1/(n^2);
    
   BZ=KZ.*KZ.*(ones(n,1)*resi.^2);
   BZ1=sum(BZ,2)-diag(BZ);
   BZ2=resi.^2*BZ1;
   BZ2=2*prod(h)*BZ2/(n^2);
 
   TfinalZ=n*sqrt(prod(h))*AZ2/sqrt(BZ2);
   
%    bootstrap
   for t=1:re
%        fprintf('t=%d\n',t);
       prob=(1+sqrt(5))/(2*sqrt(5));Z=random('binomial',1,(1-prob)*ones(1,n));
       a1=0.5*(1-sqrt(5));a2=0.5*(1+sqrt(5));
       V=(a2-a1)*Z+a1;
       err_star=resi.*V;
       y_star=beta(1)+beta(2)*x(1,:)+beta(3)*x(2,:)+beta(4)*z+beta(5)*z.^2+err_star;
%        y_star=beta(1)+beta(2)*x(1,:)+beta(3)*x(2,:)+beta(4)*z+beta(5)*z.^2+beta(6)*z(1,:).^3+err_star;
       
%        fun1=@(beta,xz)beta(1)+beta(2)*xz(1,:)+beta(3)*xz(2,:)+beta(4)*xz(3,:)+beta(5)*xz(3,:).^2;
%        fun=@(beta,xz)beta(1)+beta(2)*xz(1,:)+beta(3)*xz(2,:)+beta(4)*xz(3,:)+beta(5)*xz(3,:).^2+beta(6)*xz(3,:).^3;
       beta_star=nlinfit(xz,y_star,@myselfmodel,beta0);
%        %beta0=[1 1 1 1 1];
%        %beta_star=nlinfit(x1,y_star,@mymodel,beta0);
% %        if (isnan(beta)==1)|(isinf(beta)==1)
% %            m=m+1;
% %            continue; 
% %        end
       resi_star=y_star-(beta_star(1)+beta_star(2)*x(1,:)+beta_star(3)*x(2,:)+beta_star(4)*z+beta_star(5)*z.^2);
%        resi_star=y_star-(beta_star(1)+beta_star(2)*x(1,:)+beta_star(3)*x(2,:)+beta_star(4)*z+beta_star(5)*z.^2+beta_star(6)*z.^3);

       AZ=KZ.*(ones(n,1)*resi_star);
       AZ1=sum(AZ,2)-diag(AZ);
       AZ2=resi_star*AZ1/(n^2);
       
       BZ=KZ.*KZ.*(ones(n,1)*resi_star.^2);
       BZ1=sum(BZ,2)-diag(BZ);
       BZ2=resi_star.^2*BZ1;
       BZ2=2*prod(h)*BZ2/(n^2);
       
       TstarZZ(t)=n*sqrt(prod(h))*AZ2/sqrt(BZ2);
   
%        AZ=KZ.*(ones(n,1)*resi_star');
%        AZ1=sum(AZ,2)-diag(AZ);
%        AZ2=resi_star'*AZ1/p^2;
%     
%        BZ=KZ.*KZ.*(ones(n,1)*resi_star.^2');
%        BZ1=sum(BZ,2)-diag(BZ);
%        BZ2=resi_star.^2'*BZ1;
    
%       TstarZZ(t)=AZ2/sqrt(2*BZ2);

   end     
   TstarZ=sort(TstarZZ);
   critic=TstarZ(1,floor((19/20)*re));
   %        powersZ(l)=(TfinalZ>critic);
   ss=ss+(TfinalZ>critic);
   %    ss=ss+(TfinalZ>1.96);
   Alpha=ss/l;
   fprintf('Alpha=%f\n',Alpha);
end
toc
lambd
ss
Alpha=ss/k
end