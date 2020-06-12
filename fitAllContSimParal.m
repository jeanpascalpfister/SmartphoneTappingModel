function Data = fitAllContSimParal(sub_list,sim_list,fileIn)
% JPP 6.3.2018
% dispatch the fitting jobs on the cluster. One job corresponds to one
% subject

n = length(Data); % number of subjects
rootlog = 'log/';
rootfit = 'fit/';

% dispatch the jobs
unix('module load matlab')
for sub=sub_list    
    filetxt = [rootlog 'SUB' num2str(sub) '.txt'];
    filefit = [rootfit 'SUB' num2str(sub) '.mat'];
    str1 = 'bsub -W "04:00" -R "rusage[mem=1024]" -N ';
    str2 = ['-o ' filetxt '-J "sub=' num2str(sub) '"'];
    str3 = ' matlab -singleCompThread -nosplash -nodesktop -r ';
    str4 = ['"ClusterFitSubject(' num2str(sub) ',[' num2str(sim_list) ',' fileIn ',' filefit ']),exit"'];
   [s,w] = unix([str1 str2 str3 str4]); 
end


% collect the fitted parameters from all the jobs