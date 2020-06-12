function ClusterFitSubject(subject,serie,sim,fileIn,fileOut)
% JPP/6.3.2018 fits the priority-based model. Tuned to Euler cluster.
% mode 25.6.2018. updated version which contains model types

load(fileIn,'Data'); % load taps

dire = 'fit2';

% creates directory if it does not exist
if not(exist(dire,'dir') == 7)
    unix(['mkdir ' dire])
end

Dat = Data{subject};
[f,Model] = FitContSimSeries(Dat,serie,sim);
f.Model = Model;

save(fileOut,'f')

end
