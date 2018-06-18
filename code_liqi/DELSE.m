function [directions, lambda] = DELSE(x,w,y)

% x the observed  continuous predictor matrix
% w the observed  discrete predictor matrix
% y the observed response vector
% y=y-mean(y);
c=length(unique(w));% compute the number of different values
tildew=unique(w);
[n p1]=size(x);
P=zeros(c,1);
M=zeros(p1);
% cov_x=cov(x);
% invcov_x=inv(cov(x));
for i=1:c
    a=find(w==tildew(i));    
    x_w=x(a,:);
    y_w=y(a);
    n_w=length(y_w);
    M_w=zeros(p1);
    cov_x_w=cov(x_w);
    for j=1:n_w               %×ó³ý´úÌæÇóÄæ
        b_w=(x_w)'*(y_w<y_w(j))/n_w;
        b_w=cov_x_w\(b_w);
        b_w=b_w/(norm(b_w)+1);
        M_w=M_w+b_w*(b_w)';
    end   
    M_w=M_w/n_w;
    M=M+n_w/n*M_w;
%    UM=y_w*ones(1,n_w)-ones(n_w,1)*y_w';
%    M=inv(cov(x_w))*(x_w)'UM;

end
[directions, eig_value] = eigs(M,p1-1);
lambda=diag(eig_value);
%directions = (cov(x_orig))^(-1/2) * C;
%lambda = eig(full(mean(Lambda,3)));
%lambda = sort(lambda,'descend');




