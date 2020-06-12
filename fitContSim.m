function [theta,NegCost,lambda,thetaforce,flag,pen] = fitContSim(Data,sim)
% JPP 18.8.2017
% 9.8.2017 NegCost. Only renaming update

ITI = Data.ITI;

switch sim
    case 1 % Model 1
        disp('sim 1: theta = (a,alpha)') % rho = exp(alpha)
        lambda = [];
        thetaforce = [NaN 1 NaN 0 1]';
        theta0 = [1 1 -4 0 1E-3]'; 
        pen = 0;
    case 2 % Model 2
        disp('sim 2: theta = (a,b,alpha)')
        lambda = [];
        thetaforce = [NaN NaN NaN 0 1]';
        theta0 = Data.fit{1}.theta;
        pen = 0;
    case 3
        disp('sim 3: theta = (a,alpha,tau)')
        lambda = [];
        thetaforce = [NaN 1 NaN 1 NaN]';        
        theta0 = [Data.fit{1}.theta(1:3);1;1E-2];
        pen = 0;
    case 4
        disp('sim 4: theta = (a,b,alpha,gamma,tau)')
        lambda = [];
        thetaforce = [NaN NaN NaN NaN NaN]';
        theta0 = Data.fit{3}.theta; 
        pen = 0;
    case 5  % Model 3
        disp('sim 5: p=1E-2, theta = (a,alpha,gamma_1,...,gamma_N)  tau1 = 50, tau8 = 1000')
        n = 8;                
        x = 20^(1/(n-1));
        lambda = 1./(50*x.^(0:(n-1))); 
        thetaforce = [NaN, 1, NaN, NaN*ones(1,n)]';
        theta0 = [Data.fit{1}.theta(1:3); zeros(n,1)];
        pen = 1E-2;
    case 6  % Model 4   
        disp('sim 6: p=1E-2, theta = (a,b,alpha,gamma_1,...,gamma_N) tau1 = 50, tau8 = 1000')
        n = 8;
        x = 20^(1/(n-1));
        lambda = 1./(50*x.^(0:(n-1)));
        thetaforce = [NaN, NaN, NaN, NaN*ones(1,n)]';
        theta0 = [Data.fit{2}.theta(1:3); zeros(n,1)]; 
        pen = 1E-2;
    case 7    
        disp('sim 7: p=1E-1, theta = (a,b,alpha,gamma_1,...,gamma_N) tau1 = 50, tau8 = 1000')
        n = 8;
        x = 20^(1/(n-1));
        lambda = 1./(50*x.^(0:(n-1)));
        thetaforce = [NaN, NaN, NaN, NaN*ones(1,n)]';
        theta0 = [Data.fit{2}.theta(1:3); zeros(n,1)]; 
        pen = 1E-1;
    case 8    
        disp('sim 8: p=1E-1, theta = (a,b,alpha,gamma_1,...,gamma_N) tau1 = 50, tau12 = 1500')
        n = 12;
        x = 30^(1/(n-1));
        lambda = 1./(50*x.^(0:(n-1)));
        thetaforce = [NaN, NaN, NaN, NaN*ones(1,n)]';
        theta0 = [Data.fit{2}.theta(1:3); zeros(n,1)]; 
        pen = 1E-1;
    case 9  % Model 5  
        disp('sim 9: p=1E-2, theta = (a,b,alpha,gamma_1,...,gamma_N) tau1 = 50, tau12 = 1500')
        n = 12;
        x = 30^(1/(n-1));
        lambda = 1./(50*x.^(0:(n-1)));
        thetaforce = [NaN, NaN, NaN, NaN*ones(1,n)]';
        theta0 = [Data.fit{2}.theta(1:3); zeros(n,1)]; 
        pen = 1E-2;
    case 10    
        disp('sim 10: p=1E-3, theta = (a,b,alpha,gamma_1,...,gamma_N) tau1 = 50, tau12 = 1500')
        n = 12;
        x = 30^(1/(n-1));
        lambda = 1./(50*x.^(0:(n-1)));
        thetaforce = [NaN, NaN, NaN, NaN*ones(1,n)]';
        theta0 = [Data.fit{2}.theta(1:3); zeros(n,1)]; 
        pen = 1E-3;

    otherwise
        disp('unknown simulation')
end


[theta,NegCost,flag] = fitCont(ITI,Model);


end

