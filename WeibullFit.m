function [lambda,k,varargout] = WeibullFit(x)
% JPP 1.5.2020
%
% Extracts the parameters of the Weibull distribution
% from a list of ITI
% ex: [lambda,k,L,BIC] = WeibullFit(x)
np = 2; % number of parameters, i.e. lambda and k
N = length(x);

logxv = log(x);
logx = mean(logxv);

F = @(k)  mean(x.^k)/(mean((x.^k).*logxv)-mean(x.^k)*logx)-k;
k = fsolve(@(k) F(k),1,optimoptions('fsolve','Display','off'));
lambda = mean(x.^k)^(1/k);

L = N*(log(k/lambda) + (k-1)*(logx-log(lambda))-1);
BIC = -2*L+np*log(N);

if nargout>2
    varargout{1} = L;
end

if nargout>3
    varargout{2} = BIC;
end

