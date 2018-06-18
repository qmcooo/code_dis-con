function power=discrete_continuous_ProDimAda(n,a, lambdaw,h)
s=300;
p1=1;
p2=2;
c1=0.2*log(n)/sqrt(n);
c2=0.2*log(n)/sqrt(n);
critic1=norminv(0.99,0,1);
critic2=norminv(0.95,0,1);
critic3=norminv(0.9,0,1);
hatq1=0;
for l=1:s
    x=1*randn(n,p1); 
    w1=randn(n,p2);
    w=(w1>0);
%     w=ones(n,1);
    z=[x w];
    error=randn(n,1); % control the signal to  noise;
    
%     y=1.0*z(:,1)+z(:,2)+z(:,3)+a*cos(0.6*pi*(1.0*z(:,1)+z(:,2)+z(:,3)))+error;  %H11
%     y=sum(z,2)/sqrt(p1)+a*cos(0.6*pi*(sum(z,2)/sqrt(p1)))+error;
%     y=1.0*z(:,1)+z(:,2)+z(:,3)+w+a*exp(-(1.0*z(:,1)+z(:,2)+z(:,3)+w).^2)+error;    %H12
%     y=sum(z,2)/sqrt(p1)+a*exp(-((sum(z,2))/sqrt(p1)).^2)+error;
%     y=1.0*z(:,1)+z(:,2)+z(:,3)+z(:,1)+z(:,2)+z(:,3)+z(:,1)+z(:,2)+a*(z(:,1)+z(:,2)+z(:,3)).^2+error;     %H13 
%     y=sum(z,2)/sqrt(p1)+a*((sum(z,2)/sqrt(p1)).^2)+error;
    y=1+z(:,2)+z(:,3)+z(:,1)+z(:,1).^2+error;          %H21
    

%---------------------etimate B_1--------------------   
   [directions, eigva11] =DELSE(x,w,y);
    lambda11 = eigva11./(1+eigva11);
    lambda11 =lambda11.^2;
    d1=length(lambda11);
    lambda12=(lambda11(1:(d1-1))+c1)./(lambda11(2:d1)+c1)-1; 
    lambda13=(lambda12(2:(d1-1))+c2)./(lambda12(1:(d1-2))+c2); 
    lambda14=(lambda13<0.5);
    if sum(lambda14)<0.5
        hatq=1;
    else
        index = find(lambda14);
        hatq= index(end);
    end
    %hatq1=hatq1+hatq;
    alpha01=directions(:,1:hatq)*sign(directions(1,1));
    x2=x*alpha01; 
    MZ=zeros(n,n);
    Tx=zeros(n,n);
%     beta0=ones(1,p1+1);
    beta0=ones(1,5);
    beta=nlinfit(z,y,@myselfmodel,beta0);
%     resi=y-(sum(z*diag(beta),2));
    resi=y-(beta(1)+beta(2)*z(:,2)+beta(3)*z(:,3)+beta(4)*z(:,1)+beta(5)*z(:,1).^2);
    
    index=x*alpha01;
    x1=x*alpha01(:,1);
    %h=2.5*std(x1)*n^(-1/(4+1));
%     h=5*n^(-1/(4+1));
    Tr=(x1*ones(1,n)-ones(n,1)*x1')/h;
        MZ=exp(-Tr.^2/2)/2;         
%         MZ=15/16*(1-Tr.^2)^2;      %kernel
        KZ=MZ+0.00001;
        
%         Tw=1-abs(w*ones(1,n)-ones(n,1)*w');   %%%%check
%         Tw=lambdaw.^Tw;

%         Tw=1-((w*ones(1,n)-ones(n,1)*w')==0);    %%%%check  unordered predictor
%         Tw=lambdaw.^Tw;
        Tw=lambdaw.^abs(w*ones(1,n)-ones(n,1)*w');  %ordered predictor
        KZ=KZ.*Tw;
        KZ=KZ-diag(diag(KZ));
        A=KZ.*(ones(n,1)*resi');
        A1=sum(A,2)-diag(A);
        Vn=resi'*A1;    
        
        B=KZ.*KZ.*(ones(n,1)*resi.^2');
        B1=sum(B,2)-diag(B);
        Sn=resi.^2'*B1; 
        
        Tindex1=zeros(n,n, hatq);
    if hatq==1
        Tfinal(l)=Vn/sqrt(2*Sn);
    else
        for i=1:hatq
          Tindex1(:,:,i)=(index(:,i)*ones(1,n)-ones(n,1)*index(:,i)').^2;            
        end
        K1=1./(sum(Tindex1,3)).^0.5;
        K1=K1.*Tw;
        K1=K1-diag(diag(K1));
        An=K1.*(ones(n,1)*resi');
        A1=sum(An,2)-diag(An);
        Vn=resi'*A1; 
        
        Tfinal(l)=Vn/sqrt(2*Sn);
    end
          
     %----------p-value-------------------------  
    powers1=(Tfinal(l)>critic1);
    powers2=(Tfinal(l)>critic2);
    powers3=(Tfinal(l)>critic3);
    powers(l,:)=[powers1 powers2 powers3 hatq];
end
power=mean(powers)