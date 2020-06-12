function [tau,g] = ComputeKernel(theta,Model)
% JPP 15.2.2019

ModelType = Model.ModelType;

switch ModelType
    case 1
        delta = theta(4);
        tau = logspace(1,log10(1000),100)';
        g = (tau>delta);
    case 2
        gam = theta(4);
        lambda = theta(5);
        tau = logspace(1,log10(10*max(1./lambda)),100)';
        g = 1-exp(-tau*lambda)*gam;
    case 3
        lambda = Model.lambda;
        gam = theta(4:end);
        tau = logspace(1,log10(10*max(1./lambda)),100)';
        g = 1-exp(-tau*lambda)*gam;
end


end

