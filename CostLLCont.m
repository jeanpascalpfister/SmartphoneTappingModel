function [cost,grad] = CostLLCont(tau,theta,Model)
% JPP 31.8.2017
% mod JPP 19.6.2018. Depends on the object Model
% mod JPP 20.6.2018. factor 2 fixed in the regulariser
%
% cost function (that includes the regulariser) as well as the gradient
%
% ModelType = 1: theta = (a,b,alpha,delta)
% ModelType = 2: theta = (a,b,alpha,gamma,lambda)
% ModelType = 3: theta = (a,b,alpha,lambda1,..,lambdan)

ModelType = Model.ModelType;

% necessary?
thetaforce = Model.thetaforce;
thetaf = thetaforce;
thetaf(isnan(thetaforce)) = 0;
f = @(x) x.*isnan(thetaforce) + thetaf;
theta = f(theta);

switch ModelType
    case {1,2}        
        cost = -LLCont(tau,theta,Model);
        if nargout > 1
             grad = -gradLLCont(tau,theta,Model);
        end
    case 3
        pen = Model.pen;
        gam = theta(4:end);
        cost = -LLCont(tau,theta,Model) + pen*sum(gam.^2);
        if nargout > 1            
            grad = -gradLLCont(tau,theta,Model) + 2*pen*[zeros(3,1);theta(4:end,1)]; %20.6. changed factor 2 + simplified expression
        end
end  

end

