function [k,theta,varargout] = GammaFit(x)
% JPP 15.10.2019
% mod 1.5.2020: added varargout
%
% Extracts the parameters of the gamma distribution
% from a list of ITI
np = 2; % number of parameters, i.e. k and theta

N = length(x);
xlog = mean(x.*log(x));
m = mean(x);
logx = mean(log(x));

% initialisation
theta0 = xlog-m*logx;
k0 = m/theta0;

% find the root
F = @(k) logx - log(m/k) - psi(k);
k = fsolve(@(k) F(k),k0,optimoptions('fsolve','Display','off'));
theta = m/k;

L = N*((k-1)*logx-m/theta-k*log(theta)-log(gamma(k))); % log-likelihood
BIC = -2*L+np*log(N);

if nargout>2
    varargout{1} = L;
end

if nargout>3
    varargout{2} = BIC;
end

