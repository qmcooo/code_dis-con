function [h,lambd]=CV_main()
% diary('ex1_miss1.mat')
n=600;
k=1000;
p1=1;
p2=2;
ss=0;
re=400;

alphabet = [0 1]; prob = [0.5 0.5];
w=randsrc(n,p2,[alphabet; prob]);
x=randn(n,p1);
error=randn(n,1);  
z=[x w];
%---------------¼ÙÉè--------------------
 beta0=[1 1 1 1 1]; 
 y=beta0(1)+beta0(2)*w(:,1)+beta0(3)*w(:,2)+beta0(4)*x(:,1)+beta0(5)*x(:,1).^2+error;    % DGP0
%  y=beta0(1)+beta0(2)*w(1,:)+beta0(3)*w(2,:)+beta0(4)*x(1,:)+beta0(5)*x(1,:).^2+0.25*x(1,:).^3+error; % DGP1
%  y=beta0(1)+beta0(2)*x(1,:)+beta0(3)*x(2,:)+beta0(4)*z(1,:)+beta0(5)*z(1,:).^2+0.5*x(1,:).*z(1,:)+0.5*x(2,:).*z(1,:).^2+error; % DGP2
%  y=beta0(1)+beta0(2)*x(1,:)+beta0(3)*x(2,:)+beta0(4)*z(1,:)+beta0(5)*z(1,:).^2+(x(1,:)+x(2,:)).*sin(4*pi*z(1,:))+error;  % DGP3

%----------------------------------------
CV=0;
cv_old=10;
a0=0;
a1=[0,0];
for h=0.01:0.1:1.1
    for lam1=0.001:0.1:1.1
        for lam2=0.001:0.1:1.1
            
%         cv_old=CV;
 %      Lº¯Êý
%   lambd=b*ones(p2,1);
% lambd=0*ones(n,1);%./10;
lambd=[lam1,lam2];
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

%   Wº¯Êý
% h=a*ones(p1,1);
% h=ones(m,1).*p^(-1/(4+m));
W=zeros(n,n*p1);
Wz=zeros(n,n*p1);
for s=1:p1
    Wx=(x(:,s)*ones(1,n)-ones(n,1)*x(:,s)')/h;
    Wz(:,(s-1)*n+1:s*n)=exp(-Wx.^2/2)/h;
end
W=Wz(:,1:n);
for s=2:p1
    W=W.*Wz(:,(s-1)*n+1:s*n);
end

    KZ=L.*W;
    AZ=KZ.*(ones(n,1)*y');
    AZ1=sum(AZ,2)-diag(AZ);
    AZ2=sum(KZ,2)-diag(KZ);
    CV=sum((y-AZ1./AZ2).^2)/n;
    if CV<=cv_old
        cv_old=CV
        a0=h
        a1=lambd
    end
        end
    end
end
h=a0
lambd=a1

end