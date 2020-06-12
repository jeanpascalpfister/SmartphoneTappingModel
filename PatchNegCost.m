function Data = PatchNegCost(Data,sim_list)
% JPP 9.10.2017
% Patch: replaces the LL field by NegCost and computes the true LL

nsub = length(Data);
m = length(sim_list);

for k=1:nsub
    for s=1:m
        sim = sim_list(s);
        disp(['sub = ' num2str(k) ' , sim =  ' num2str(sim)])
        tau = Data{k}.ITI;
        n = length(tau);
        Df = Data{k}.fit{sim};
        theta = Df.theta;
        lambda = Df.lambda;
        if isfield(Df,'pen')
            pen = Df.pen;
        else
            pen = 0;
        end
        
        LLn = LLCont(tau,theta,lambda); % normalised log-likelihood
        LL = LLn*n;
        NegCost = -LL + pen*sum(theta(4:end).^2);                
        Npar = length(find(isnan(Df.thetaforce)));       
        BIC = Npar*log(n) + 2*NegCost; 
        
        
        Data{k}.fit{sim}.NegCost = NegCost;
        Data{k}.fit{sim}.LL = LL;
        Data{k}.fit{sim}.LLn = LLn;
        Data{k}.fit{sim}.Npar = Npar;
        Data{k}.fit{sim}.BIC = BIC;
    end
end

end

