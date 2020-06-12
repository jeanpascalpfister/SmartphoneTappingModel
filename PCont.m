function P = PCont(tau,theta,Model)
% JPP 8.8.2017
% mod JPP 19.6.2018. Includes ModelTypes
%
% % computes the ITI distribution for the continuous tapping model with
% parameter theta
%
% ModelType = 1: theta = (a,b,alpha,delta)
% ModelType = 2: theta = (a,b,alpha,gamma,lambda)
% ModelType = 3: theta = (a,b,alpha,lambda1,..,lambdan)
%
% P = \rho_\tau\int_0^1 x exp(-r(\tau)x)p(x)dx

a = theta(1);
b = theta(2);

rho = rhotau(tau,theta,Model);
r = rtau(tau,theta,Model);

lastwarn('');
Q = integral(@(x) (exp(-x*r)-exp(-r))*betapdf(x,a+1,b),0,1,'ArrayValued',true);    

if strncmpi(lastwarn,'Infinite',length('Infinite'))
    error('JPP stop')
end

P = rho.*(Q+exp(-r))*beta(a+1,b)/beta(a,b);

end

 