function [mu,sig2,varargout] = LogNormalFit(x)
% JPP 1.5.2020
%
% Extracts the parameters of the log-normal distribution
% from a list of ITI
np = 2; % number of parameters, i.e. mu and sigma^2

N = length(x);

mu = mean(log(x));
sig2 = mean((log(x)-mu).^2);

L = -N*(mu + log(2*pi*sig2)/2 + 1/2);

BIC = -2*L+np*log(N);

if nargout>2
    varargout{1} = L;
end

if nargout>3
    varargout{2} = BIC;
end

