function [theta,NegCost,flag] = fitCont(ITI,Model)
% JPP 8.8.2017
% mod JPP 19.6.2018. Includes ModelTypes. 
% all Model options are in Model
%
% fits the continuous model where
%
% ModelType = 1: theta = (a,b,alpha,delta)
% ModelType = 2: theta = (a,b,alpha,gamma,lambda)
% ModelType = 3: theta = (a,b,alpha,lambda1,..,lambdan)
%
% Ex: [theta,LL] = fitCont(Data{1}.ITI,Model)
%
% previously. Compatibility TBC
% thetaforce = varargin{1};
% theta0 = varargin{2};
% pen = varargin{3}

ModelType = Model.ModelType;
thetaforce = Model.thetaforce; % note that 
theta0 = Model.theta0;
FitAlgo = Model.FitAlgo;
        
switch ModelType
    case 1
        delta0 = theta0(end); % should not be optimised
        thetalb = [0.05,  0.05,-20   ,delta0]'; % TBC if KO that lb = ub for delta
        thetaup = [10,      10,  0   ,delta0]';                   
        A = [];
        B = [];
    case 2
      
        thetalb = [0.05,  0.05,-20   ,-100  ,1E-8]'; 
        thetaup = [10,      10,  0   ,1     ,1   ]';                   
        A = [];
        B = [];

    case 3
        lambda = Model.lambda;
        n = length(lambda);       
        thetalb = [0.05,  0.05, -20   , -200*ones(1,n)]'; 
        thetaup = [10,      10,   0,    200*ones(1,n)]';  
        
        %N = 70;
        %t = round([0:5:40 45:65 logspace(log10(66),log10(3/lambda(end)),N) ])'; 
        
        t = round([0:1:99 logspace(log10(100),log10(3/lambda(end)),100)])'; 
        
        N1 = length(t);
        A = [zeros(N1,3) exp(-t*lambda)];
        B = ones(N1,1);
end


thetaf = thetaforce;
thetaf(isnan(thetaforce)) = 0;
f = @(x) x.*isnan(thetaforce) + thetaf;

% note that fmincon of matlab is not working well for this cost function 
% (sometimes the cost function is higher at the end of the  minimization!!!). 
% This is why I wrote my own optimisation function.

opt = optimoptions('fmincon','Display','iter','TolFun',1E-6,'TolX',1E-5,'MaxIter',1000);    

if strcmp(FitAlgo,'JPP')
    [theta,Cost,flag] = fminconJPP(@(x) CostLLCont(ITI,x,Model),theta0,A,B,[],[],thetalb,thetaup,[],opt);
elseif strcmp(FitAlgo,'matlab')
    [theta,Cost,flag] = fmincon(@(x) CostLLCont(ITI,x,Model),theta0,A,B,[],[],thetalb,thetaup,[],opt);
else
    error('unknown fitting algorithm')
end

theta = f(theta); 
NegCost = -Cost;

switch ModelType
    case 1
         disp(['ML fitting: C = ' num2str(Cost) ', a = ' num2str(theta(1)) ', b = ' num2str(theta(2)) ', alph = ' num2str(theta(3)) ', delta = ' num2str(theta(4))])
    case 2
         disp(['ML fitting: C = ' num2str(Cost) ', a = ' num2str(theta(1)) ', b = ' num2str(theta(2)) ', alph = ' num2str(theta(3)) ', gamma = ' num2str(theta(4)) ', lambda = ' num2str(theta(5))])
    case 3
         disp(['ML fitting: C = ' num2str(Cost) ', a = ' num2str(theta(1)) ', b = ' num2str(theta(2)) ', alph = ' num2str(theta(3)) ', gamma = ' num2str(theta(4:end)')])
end



