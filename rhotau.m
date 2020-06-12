function rh = rhotau(tau,theta,Model)
% JPP 8.8.2017
% mod JPP 15.8.2017. added varargin (i.e. multiple time constants)
% mod JPP 19.6.2018.  Depends on the object Model
%
% ModelType = 1: theta = (a,b,alpha,delta)
% ModelType = 2: theta = (a,b,alpha,gamma,lambda)
% ModelType = 3: theta = (a,b,alpha,lambda1,..,lambdan)

ModelType = Model.ModelType;
rho = exp(theta(3));

switch ModelType
    case 1
        delta = theta(4);
        rh = rho*(tau>delta);
    case 2
        gam = theta(4);
        lambda = theta(5);
        rh = rho*(1-exp(-tau*lambda)*gam);
    case 3
        lambda = Model.lambda;
        gam = theta(4:end);
        rh = rho*(1-exp(-tau*lambda)*gam);
end


if sum(rh < 0)>0
    disp('JPP: rectifying the rates')
    disp(['min(rho) = ' num2str(min(rh))])
    %error('error JPP: negative rate')
    rh = rh.*(rh>0);
end

