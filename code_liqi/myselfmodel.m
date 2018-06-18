function yhat=myselfmodel(beta,z)

%   yhat=sum(z*diag(beta),2);
  yhat=beta(1)+beta(2)*z(:,2)+beta(3)*z(:,3)+beta(4)*z(:,1)+beta(5)*z(:,1).^2;    %liqi
end