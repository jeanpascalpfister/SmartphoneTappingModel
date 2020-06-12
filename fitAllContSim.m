function Data = fitAllContSim(DataIn,varargin)
% JPP 12.1.2018
% Data = fitAllContSim(DataIn,sims,varargin)


if nargin>1
    file = varargin{1};    
end

n = length(DataIn);
Data = DataIn;
disp('ML Cont Sim fitting')
serie_list = 1:6;
nsim = [11 11 1 1 12 12]; % depends on the fact that delta_list contains 11 elements

for sub=1:n
    disp(['subject = ' num2str(sub)])        
    for ss=1:6
        serie = serie_list(ss);
        for sim=1:nsim(serie)        
            [f,Model] = FitContSimSeries(Data{sub},serie,sim);
            
             if nargin>2
                disp(['save serie ' num2str(serie) ', simulation ' num2str(sim) ', subject ' num2str(sub)])          
                save(file,'Data')
             end        
           
            Data{sub}.SimSerie{serie}.fit{sim} = f;
            Data{sub}.SimSerie{serie}.fit{sim}.Model = Model;
            
        end
    end
end


end

