function [taus,exitflag] = taustar(theta,Model)
% JPP 5.4.2019
%
% computes the interval taus at which the r(taus)/r(infinity) = 0.5

if Model.ModelType == 1
   error('JPP: taustar can not be computed for Model.Type 1') 
end

rho = exp(theta(3));

%tau0 = 100;%
options = optimoptions('fsolve','Display','off');

rrel = @(tau) rhotau(tau,theta,Model)/rho-0.5;

tau0list = [10 100 500 1000];
for k=1:length(tau0list)
    tau0 = tau0list(k);
    [tauslist(k),fval,exitflaglist(k)] = fsolve(rrel,tau0,options);
end

ind = find(exitflaglist==1);
if length(ind)>=1
    taus = tauslist(ind(end));
    exitflag = exitflaglist(ind(end));    
else
    %taus = NaN;
    taus = tauslist(1);
    exitflag = exitflaglist(1);    
end

%disp(num2str(exitflag))
%if exitflag ~= 1
%    disp(['******* ' num2str(exitflag)])
%else
    %disp(num2str(exitflag))
%end

end

