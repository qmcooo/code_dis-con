function lead_discrete_continuous_ProDimAda
diary('ex1_miss1.mat')

'---------lead_discrete_continuous_ProDimAda-----'

num=[100];
% b=[8:12]*0.15;
% a=[0 0.2 0.4 0.6 0.8 1.0];
% a=[0 0.03 0.05 0.08 0.1 0.2 0.35 0.5 1 2];
a=0;
resu=zeros(length(a),4);
%lambdaw=0.05;
% lambdaw=5*num^(-2/(4+1));
% h=1.5*num^(-1/(4+1));
%  lambdaw=0.1;
%  h=0.1;
% 0.0001:0.1:1.1
a11=[];
a22=[];
tt=[];
for a1=0.0001:0.1:1.1
    for a2=0.0001:0.1:1.1
        lambdaw=a1;
        h=a2;
        a1
        a2
        

for i=1:length(a)
% for j=1:length(num)
     power=discrete_continuous_ProDimAda(num, a(i),lambdaw,h);
     resu(i,:)=power
% end
end
    if abs(resu(2)-0.05)<=1e-2
        a11=[a11 a1];
        a22=[a22 a2];
        tt=[tt resu(2)];
%         flag=1;
%         break;
%     end
    end
%     if flag==1
%         break
    end
    end
dlmwrite('H21_a11.txt', a11);
dlmwrite('H21_a22.txt', a22);
dlmwrite('tt.txt', tt);
a11
a22

all_resu=resu;
dlmwrite('ex1_miss1.txt',all_resu,'delimiter',' ','newline','pc') 