function [f,Model] = FitContSimSeries(Data,serie,sim)
% JPP 25.6.2018
% 
% Fits the model for a given simulation within serie


ITI = Data.ITI;

% typical initial conditions
a0 = 1;
b0 = 1;
alpha0 = -4;
gamma0 = 0;
lambda0 = 1E-3;

delta_list = 0:5:50; % list of hard refractory delay (case 1 and 2)

switch serie
    case 1 
        % theta = (a,b=1,alpha,delta)
        Model.ModelType = 1;
        b = 1;        
        delta = delta_list(sim);        
        Model.thetaforce =  [NaN b NaN delta]';
        Model.theta0 = [a0 b alpha0 delta]';
        Model.pen = 0;
        Model.FitAlgo = 'matlab';
    case 2
        % theta = (a,b,alpha,delta)        
        Model.ModelType = 1;
        delta = delta_list(sim);                
        Model.thetaforce =  [NaN NaN NaN delta]';
        Model.theta0 = [a0 b0 alpha0 delta]';
        Model.pen = 0;
        Model.FitAlgo = 'matlab';
    case 3
        % theta = (a,b=1,alpha,gamma,lambda)
        Model.ModelType = 2;
        b = 1;
        Model.thetaforce = [NaN b NaN NaN NaN]';
        Model.theta0 = [a0 b alpha0 gamma0 lambda0]';
        Model.pen = 0; 
        Model.FitAlgo = 'matlab';
    case 4
        % theta = (a,b,alpha,gamma,lambda)
        Model.ModelType = 2;       
        Model.thetaforce = [NaN NaN NaN NaN NaN]';
        Model.theta0 = [a0 b0 alpha0 gamma0 lambda0]';
        Model.pen = 0;                 
        Model.FitAlgo = 'matlab';
    case 5
        % theta = (a,b=1,alpha,gamma1,...,gamman)           
        Model.ModelType = 3;       
        b = 1;
        n = sim;
        x = 20^(1/(n-1));        
        Model.lambda = 1./(50*x.^(0:(n-1)));
        Model.thetaforce = [NaN b NaN NaN*ones(1,n)]';
        Model.theta0 = [a0 b alpha0 zeros(1,n)]';
        Model.pen = 1E3; % changed from 1E-2                 
        %Model.FitAlgo = 'JPP';
        Model.FitAlgo = 'matlab';
    case 6
        % theta = (a,b,alpha,gamma1,...,gamman)           
        Model.ModelType = 3;               
        n = sim;
        x = 20^(1/(n-1));        
        Model.lambda = 1./(50*x.^(0:(n-1)));
        Model.thetaforce = [NaN NaN NaN NaN*ones(1,n)]';
        Model.theta0 = [a0 b0 alpha0 zeros(1,n)]';
        Model.pen = 1E3; % changed from 1E-2                 
        %Model.FitAlgo = 'JPP';
        Model.FitAlgo = 'matlab';
end
     
%Model.FitAlgo = 'JPP';
[theta,NegCost,flag] = fitCont(ITI,Model);

n = length(ITI);
f.theta = theta;
f.NegCost = NegCost;
f.flag = flag;
f.LL = LLCont(ITI,theta,Model);
f.LLn = f.LL/n;
f.Npar = length(find(isnan(Model.thetaforce)));   
f.BIC = f.Npar*log(n) - 2*f.NegCost; % + => - changed on 23.8.2018


%        LLn = LLCont(tau,theta,lambda); % normalised log-likelihood
%        LL = LLn*n;
%        NegCost = -LL + pen*sum(theta(4:end).^2);                
%        Npar = length(find(isnan(Df.thetaforce)));       
%        BIC = Npar*log(n) + 2*NegCost; 
        

        
end


