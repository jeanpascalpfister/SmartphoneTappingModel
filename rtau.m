function r = rtau(tau,theta,Model)
% JPP 8.8.2017
% mod 19.6.2018. Depends on the object Model
% 
% computes the integral of the tapping rate from 0 to tau
%
% ModelType = 1: theta = (a,b,alpha,delta)
% ModelType = 2: theta = (a,b,alpha,gamma,lambda)
% ModelType = 3: theta = (a,b,alpha,lambda1,..,lambdan)

ModelType = Model.ModelType;

rho = exp(theta(3));

switch ModelType
    case 1
        delta = theta(4);
        r = rho*(tau-delta).*(tau>delta);
    case 2
        gam = theta(4);
        lambda = theta(5);
        r = rho*(tau-(1-exp(-tau*lambda))*(gam./lambda'));
    case 3
        lambda = Model.lambda;
        gam = theta(4:end);
        r = rho*(tau-(1-exp(-tau*lambda))*(gam./lambda'));
end
         
if sum(r < 0)>0
    
    disp('JPP: rectifying the rates')
    error('error JPP: negative integral of rate')
    r = r.*(r>0);
    
end


end

