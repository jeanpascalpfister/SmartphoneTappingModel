function dtheta = gradLLCont(tau,theta,Model)
% JPP 30.8.2017
% mod 19.6.2018. Includes ModelTypes
% mod 23.8.2018 mean->sum
% 
% computes the gradient of the LL
%
% ModelType = 1: theta = (a,b,alpha,delta)
% ModelType = 2: theta = (a,b,alpha,gamma,lambda)
% ModelType = 3: theta = (a,b,alpha,lambda1,..,lambdan)

ModelType = Model.ModelType;
thetaforce = Model.thetaforce;

a = theta(1);
b = theta(2);
rho0 = exp(theta(3));

% lambda,thetaforce,ModelType

rho = rhotau(tau,theta,Model);
r = rtau(tau,theta,Model);

E1 = exp(-r);

% <xE(x)>
m = a/(a+b);
Q = integral(@(x) (exp(-x*r)-E1)*betapdf(x,a+1,b),0,1,'ArrayValued',true);    
xE = m*(Q+E1);

% <log(x)>
Q1 = integral(@(x) x^(a-1)*((1-x)^(b-1)-(1-x))*log(x),0,1,'ArrayValued',true);    
c = 1/(a+1)^2-1/a^2;
logx = (Q1+c)/beta(a,b);

% <log(1-x)>_(a,b)
d = 1/(b+1)^2-1/b^2;
Q2 = integral(@(x) (x^(a-1)-x)*(1-x)^(b-1)*log(1-x),0,1,'ArrayValued',true);
log1x = (d+Q2)/beta(a,b);

% <log(1-x)>_(a+1,b)
Q2a1 = integral(@(x) (x^(a)-x)*(1-x)^(b-1)*log(1-x),0,1,'ArrayValued',true);
log1xa1 = (d+Q2a1)/beta(a+1,b);


%<xE(x)log(x)>
Q3 = integral(@(x) exp(-x*r)*betapdf(x,a+1,b)*log(x),0,1,'ArrayValued',true);    
xElogx = Q3*m;

%<xE(x)log(1-x)>
Q4 = integral(@(x) (exp(-x*r)-E1)*log(1-x)*betapdf(x,a+1,b),0,1,'ArrayValued',true);    
xElog1x = m*(Q4+E1*log1xa1);

%<x^2E(x)>
Q5 = integral(@(x) (exp(-x*r)-E1)*betapdf(x,a+2,b),0,1,'ArrayValued',true);    
x2E = (Q5+E1)*(beta(a+2,b)/beta(a,b));


% gradient
da   = sum(xElogx./xE-logx);        % dL/da
db   = sum(xElog1x./xE-log1x);      % dL/db
%drho = mean((1-(x2E./xE).*r)./rho); % dL/drho % !!! rho or rho0?
dalph = sum(1-(x2E./xE).*r);        %dL/dalpha % change of param

switch ModelType
    case 1
        %delta = theta(5); 
        ddelta = 0;% since the LL is not differentiable w.r.t. delta, no gradient computation here
        dtheta = [da;db;dalph;ddelta];                
    case 2
        gam = theta(4);
        lambda = theta(5);
        
        n = 1;
        %B: basis functions evaluated at all ITI's
        B = exp(-tau*lambda);
        M1 = -B./repmat(rho,1,n);
        M2 = -((x2E./xE)*(1./lambda)).*(1-B);
        
        dlam = rho0*gam*sum((B.*tau)./rho + (x2E./xE).*((B-1)/lambda^2+B.*tau/lambda));                    
        dgam = rho0*sum(M1-M2,1)';          % dL/dgamma
        dtheta = [da;db;dalph;dgam;dlam];        
    case 3
        %gam = theta(4:end); % required?
        lambda = Model.lambda;
        n = length(lambda);
        B = exp(-tau*lambda);
        M1 = -B./repmat(rho,1,n);
        M2 = -((x2E./xE)*(1./lambda)).*(1-B);
        
        dgam = rho0*sum(M1-M2,1)';          % dL/dgamma
        dtheta = [da;db;dalph;dgam];       
end

dtheta = dtheta.*isnan(thetaforce);

end

