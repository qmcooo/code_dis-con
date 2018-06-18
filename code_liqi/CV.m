function [h0,lambd0]=CV(x,y,z,n,p,q)
CV11=n*n;

for a=0.01:0.05:2.01
 %   for b=0.05:0.1:0.55
b=0;
h1=a*n^(-1/(4+q));
lambd1=b*n^(-2/(4+q));
        %      Lº¯Êý
lambd=lambd1*ones(p,1);
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
h=ones(q,1).*h1;
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

KZ=L.*W;
AZ=((sum(KZ,2)-diag(KZ))./n).^-1;
AZ1=KZ*y'-diag(KZ).*y';
AZ2=(AZ.*AZ1)./n;

CV=sum((y'-AZ2).^2)/n;
if CV<CV11
    a1=a;
    b1=b;
    h0=h1;
    lambd0=lambd1;
    CV11=CV;
end
end
%b1
% h0
% lambd0
% CV1
end