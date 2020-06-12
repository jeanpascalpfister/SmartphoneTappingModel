function LL = LLCont(tau,theta,Model)
% JPP 8.8.2017
%mod JPP 19.6.2018.  Depends on the object Model
% mod JPP 23.6.2018 removed the normalisation
%
% computes the (unnormalised) log-likelihood for the continuous tapping model with
% parameter theta
% theta = (a,b,rho,gamma,lambda)
%
% P = \rho_\tau\int_0^1 x exp(-r(\tau)x)p(x)dx

%n = length(tau);
LL = sum(log(PCont(tau,theta,Model))); %/n;


end