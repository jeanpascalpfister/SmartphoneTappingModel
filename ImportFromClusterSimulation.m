function Data = ImportFromClusterSimulation(fileIn)
% JPP 28.8.2017

load(fileIn,'Data');

nsim = [11 11 1 1 50 50]; % could be improved by checking the number of simulations per serie

for serie=1:6
    
    for sim = 1:nsim(serie)
        
        for sub=1:length(Data)
        
            file = ['fit2/SER' num2str(serie) '_SUB' num2str(sub) '_SIM' num2str(sim) '.mat'];                
    
            if exist(file,'file') == 2 
                load(file,'f'); 
                Data{sub}.SimSerie{serie}.fit{sim} = f;       
            else
                %disp(['JPP:  subject ' num2str(sub) ' not downloaded'])
                disp(['JPP: file ' file ' not found'])
            end
        end
    end
end
    