function GenerateFiguresRevision1
% JPP 1.3.2017
%
% generates all the figures for the following paper:
% Pfister, J. P., & Ghosh, A.  To touch or not to touch: a multi-scale priority model for smartphone actions.

% file names
DataFile = 'DataSet.mat';           % tapping data
FitFile = 'DataAndModelFit.mat'; % tapping data + model parameter fits
root = 'figures/';                  % folder in which figures will be saved
regen = 0;                          % 1: regenerates the fitting, 0: load fitted results from file
parallelize = 1;                    % 1: if the simulation can be run on a cluster, 0: else
postprocess = 0;

% fitting options
sub =10;  % subject number
SimSerie = 6; fit = 21;           % relative refractoriness best fit
SimSerieHard = 2; fitHard = 11;   % hard refractoriness best fit

% figures options
fsize = 40; % default Font Size

if not(exist(root,'dir') == 7)
    unix(['mkdir ' root])
end

% unix('bash cropall.sh') % crops all the pdfs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        MODEL FITTING                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(FitFile) == 2 && regen == 0
    disp('** load smartphone touching models  **')
    load(FitFile); % loads Data
else    
    disp('** ML Fitting of the continuous smartphone touching models  **')    
    if exist(FitFile) == 2    
        load(FitFile); 
    else
        load(DataFile); % loads Data
    end
    if parallelize == 1
        disp('1. Copy all the files on the cluster')
        disp('2. Run the shell script: RunCluster.sh')
        disp('3. Wait until the simulation is finished')
        disp('4. copy the fit2/ directory from the cluster')
        struser = input('Did you perform those 4 steps? [yes/no]','s');
        if strcmp(struser,'yes')
            %Import simulation results from the cluster
            Data = ImportFromClusterSimulation(DataFile);
        else            
            disp('Simulations have not been run on the cluster')
            return; 
        end
    else
        Data = fitAllContSim(Data,FitFile);
    end
    %Data = PatchNegCost(Data,sim_list); % compute BIC
    save(FitFile,'Data');
end

if postprocess == 1

    Nf = length(Data{1}.SimSerie{5}.fit);
    nexit = 0;
    for k=1:length(Data)
        for ss=5:6
            for f=1:Nf
                Model = Data{k}.SimSerie{ss}.fit{f}.Model;
                theta = Data{k}.SimSerie{ss}.fit{f}.theta;                         
                [taus, exitflag] = taustar(theta,Model); 
                Data{k}.SimSerie{ss}.fit{f}.taus = taus;
                Data{k}.SimSerie{ss}.fit{f}.tausexitflag = exitflag;

                if exitflag ~= 1
                    disp(['exitflag = ' num2str(exitflag) ', sub = ' num2str(k) ', ss = ' num2str(ss) ', f = ' num2str(f)])
                    nexit = nexit + 1;
                end
            end       
        end
    end
    disp(['n exit = ' num2str(nexit)])
    save(FitFile,'Data');
end

Model = Data{sub}.SimSerie{SimSerie}.fit{fit}.Model;
ModelHard = Data{sub}.SimSerie{SimSerieHard}.fit{fitHard}.Model;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              FIGURES                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set FigXX = 1 to display and save the figure

% Figure 1 "Smartphone touch data"
Fig1a = 0;
Fig1b = 0;

% Figure 2 "Properties of the smartphone touching model"
Fig2a = 0;
Fig2b = 0;

% Figure 3 "Fitting results for the models with hard refractoriness"
Fig3a = 0;
Fig3b = 0;
Fig3c = 0;
Fig3d = 0;

% Figure 4 "Fitting results of model M6"
Fig4a = 0;
Fig4b = 0;
Fig4c = 0;
Fig4d = 0;
Fig4e = 0;
Fig4f = 0;
Fig4g = 0;
Fig4h = 0;

% Figure 5 "Model Comparison"
Fig5a = 1;
Fig5b = 1;
Fig5c = 1;

% Figure 6 "The scale-free exponent a is inversely correlated with the effective ref. time cst"
Fig6a = 0;
Fig6b = 1;
Fig6binset = 0;
Fig6c = 0;
Fig6d = 0;
Fig6e = 0;
Fig6f = 0;


% % OLD FIGURE NAMES %
% 
% % Figure 0: data
% Fig0A = 0; %paper fig1a
% Fig0E = 0; %paper fig1b
% 
% % Figure 1: model
% Fig1G = 0; %paper fig2b
% Fig1H = 0; %paper fig2a
% 
% % Figure 2: results
% Fig2A = 0; %paper fig4a1 (when SimSerie=6,fit=21) and fig3c (when SimSerie=2,fit=11)
% Fig2B = 1; %paper fig4b1 
% Fig2C = 0; %paper fig4c1
% Fig2D = 0; %paper fig4a2 (when SimSerie=6,fit=21) and fig3d (when SimSerie=2,fit=11)
% Fig2E = 0; %paper fig4b2
% Fig2F = 0; %paper fig4c2
% Fig2G = 0; %paper fig4d
% Fig2H = 0; %paper fig4e
% 
% % Figure 3: model comparison
% Fig4A = 0; % paper fig3a.  Log-likelihood as a function of Delta for one subject
% Fig4B = 0; % paper fig3b  log-likelihood as a function of Delta across all subjects
% Fig4Db = 0; % paper fig5a.  % model comparison for the whole population across all models BICtot
% Fig4Dc = 0; % paper fig5a. new  % figure for revised version of the paper.  
% Fig4F = 0; % paper fig5b
% Fig4G = 0; % paper fig5c
% Fig4I = 0; % paper fig 6a-f


prefix = '';

filenameFig1a = [prefix 'Fig1a'];
filenameFig1b = [prefix 'Fig1b'];

filenameFig2a = [prefix 'Fig2a'];
filenameFig2b = [prefix 'Fig2b'];

filenameFig3a = [prefix 'Fig3a'];
filenameFig3b = [prefix 'Fig3b'];
filenameFig3c = [prefix 'Fig3c'];
filenameFig3d = [prefix 'Fig3d'];

filenameFig4a = [prefix 'Fig4a'];
filenameFig4b = [prefix 'Fig4b'];
filenameFig4c = [prefix 'Fig4c'];
filenameFig4d = [prefix 'Fig4d'];
filenameFig4e = [prefix 'Fig4e'];
filenameFig4f = [prefix 'Fig4f'];
filenameFig4g = [prefix 'Fig4g'];
filenameFig4h = [prefix 'Fig4h'];

filenameFig5a = [prefix 'Fig5a'];
filenameFig5b = [prefix 'Fig5b'];
filenameFig5c = [prefix 'Fig5c'];

filenameFig6a = [prefix 'Fig6a'];
filenameFig6b = [prefix 'Fig6b'];
filenameFig6c = [prefix 'Fig6c'];
filenameFig6d = [prefix 'Fig6d'];
filenameFig6e = [prefix 'Fig6e'];
filenameFig6f = [prefix 'Fig6f'];


% filenameFig0B = [prefix 'Fig0B'];
% filenameFig0C = [prefix 'Fig0C'];
% filenameFig0D = [prefix 'Fig0D'];
% filenameFig0E = [prefix 'Fig0E'];
% 
% filenameFig1G = [prefix 'Fig1G'];
% filenameFig1H = [prefix 'Fig1H'];
% 
% filenameFig2A = [prefix 'Fig2A'];
% filenameFig2B = [prefix 'Fig2B'];
% filenameFig2C = [prefix 'Fig2C'];
% filenameFig2D = [prefix 'Fig2D'];
% filenameFig2E = [prefix 'Fig2E'];
% filenameFig2F = [prefix 'Fig2F'];
% filenameFig2G = [prefix 'Fig2G'];
% filenameFig2H = [prefix 'Fig2H'];
% 
% filenameFig4A = [prefix 'Fig4A'];
% filenameFig4B = [prefix 'Fig4B'];
% filenameFig4Db = [prefix 'Fig4Db'];
% filenameFig4Db1 = [prefix 'Fig4Db1'];
% filenameFig4Dc = [prefix 'Fig4Dc'];
% %filenameFig4E = [prefix 'Fig4E'];
% filenameFig4F = [prefix 'Fig4F'];
% filenameFig4G = [prefix 'Fig4G'];
% filenameFig4I = [prefix 'Fig4H'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   FIGURE 1: Smartphone touch data   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Fig1a == 1

    [f1,f2,f3,f4] = PlotNiceTaps(Data);
    
    print(f1,'-dpdf',[root filenameFig1a '1.pdf'])
    savefig(f1,[root filenameFig1a '1.fig'])
    
    print(f2,'-dpdf',[root filenameFig1a '2.pdf'])
    savefig(f2,[root filenameFig1a '2.fig'])

    print(f3,'-dpdf',[root filenameFig1a '3.pdf'])
    savefig(f3,[root filenameFig1a '3.fig'])

    print(f4,'-dpdf',[root filenameFig1a '4.pdf'])
    savefig(f4,[root filenameFig1a '4.fig'])

end

if Fig1b == 1
    f5 = setfigure(5);
    clf;
    
    ITI = Data{sub}.ITI;
    PlotITI(ITI); % Data    
    setaxis2(fsize)
    xlabel('ITI $\tau$ [ms]','interpreter','latex')
    ylabel('$p(\tau)$','interpreter','latex')
    xlim([50 max(ITI)])
    set(gca,'Xtick',[1E2 1E4 1E6 1E8])
    
    %PlotITIContDataFit(Data,sub,theta,Model,0)
    print(f5,'-dpdf',[root filenameFig1b '.pdf'])
    savefig(f5,[root filenameFig1b '.fig'])
end  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIGURE 2: Model Properties          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Fig2a ==1
    listtaur = [100 300 1000];
    col = [0 0 0; 0.5 0 0; 1 0 0];
    style = {'-','-.',':'};
    tau = logspace(log10(10),8,100)';
    rho = 1E-3; a = 1; b = 1; gam = 1;
    alph = log(rho);
    
       
    % Model definition
    Mod.ModelType = 2;
    Mod.thetaforce = [NaN NaN NaN NaN NaN];
    Mod.theta0 = [1 1 1 1 1];
    Mod.pen = 0;
    Mod.FitAlgo = 'matlab';
    
    f41 = setfigure(21);
    clf;
    for k=1:length(listtaur)
        lam = 1/listtaur(k);
        theta = [a,b,alph,gam,lam];
        p(:,k) = PCont(tau,theta,Mod);
        h(k) = loglog(tau,p(:,k),'color',col(k,:),'linewidth',2,'linestyle',style{k});
        LEG{k} = ['$\tau_r$ = ' num2str(1/lam) ' ms'];
        if k==1
            hold on
        end                
    end
    
    xlim([1E1 1E8])
    ylim([1E-12 1E-2])
    set(gca,'Xtick',[1E2 1E4 1E6 1E8],'Ytick',[1E-12 1E-10 1E-8 1E-6 1E-4 1E-2])
    legend(h,LEG,'box','off','interpreter','latex','position',[0.67 0.7 0.2 0.2]);
    setaxis2(fsize)    
    
    xlabel('ITI $\tau$ [ms]','interpreter','latex')
    ylabel('$p(\tau)$','interpreter','latex')
    
    print(f41,'-dpdf',[root filenameFig2a '.pdf'])
    savefig(f41,[root filenameFig2a '.fig'])
    
end


if Fig2b ==1
    lista = [1 1.5 2];
    col = [0 0 0; 0.5 0 0; 1 0 0];
	style = {'-','-.',':'};
    tau = logspace(log10(10),8,100)';
    rho = 1E-3; b = 1; gam = 1;lam = 1/100;
    alph = log(rho);    
    
    % Model definition
    Mod.ModelType = 2;
    Mod.thetaforce = [NaN NaN NaN NaN NaN];
    Mod.theta0 = [1 1 1 1 1];
    Mod.pen = 0;
    Mod.FitAlgo = 'matlab';
    
    f43 = setfigure(22);
    clf;
    for k=1:length(lista)
        a = lista(k);
        theta = [a,b,alph,gam,lam];
        p(:,k) = PCont(tau,theta,Mod);
        h(k) = loglog(tau,p(:,k),'color',col(k,:),'linewidth',2,'linestyle',style{k});
        LEG{k} = ['$a$ = ' num2str(a)];
        if k==1
            hold on
        end                
    end
    
    xlim([1E1 1E8])
    ylim([1E-12 1E-2])
    set(gca,'Xtick',[1E2 1E4 1E6 1E8],'Ytick',[1E-12 1E-10 1E-8 1E-6 1E-4 1E-2])
    legend(h,LEG,'box','off','interpreter','latex');    
    setaxis2(fsize)    
    
    xlabel('ITI $\tau$ [ms]','interpreter','latex')
    ylabel('$p(\tau)$','interpreter','latex')    
    
    print(f43,'-dpdf',[root filenameFig2b '.pdf'])
    savefig(f43,[root filenameFig2b '.fig'])
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIGURE 3: Fitting results hard ref  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LL as a function of \Delta
if Fig3a == 1
   
    % for a single subject
    for k=1:11        
        LL1(k) = Data{sub}.SimSerie{1}.fit{k}.LL;
        LL2(k) = Data{sub}.SimSerie{2}.fit{k}.LL;       
    end
   
    f31 = figure(31);
    clf;
    delta = 0:5:50;   
    
    l(1) = plot(delta,LL2,'r-','linewidth',2);
    hold on
    l(2) = plot(delta,LL1,'k-.','linewidth',2);    
    
    
    setaxis2(fsize)
    set(gca,'Xtick',[0 20 40]);
    legend(l,{'$b\neq 1$','$b=1$'},'box','off','interpreter','latex','Location','NorthWest')
    
    
    xlabel('refractory time $\Delta$ [ms]','interpreter','latex')
    ylabel('log-likelihood $L$','interpreter','latex')
    xlim([0 50])
    
    print(f31,'-dpdf',[root filenameFig3a '.pdf'])
    savefig(f31,[root filenameFig3a '.fig'])
    
end

% LL as a function of \Delta pooled across subjects
if Fig3b == 1
    
    for subj = 1:length(Data)
        ntaps = length(Data{subj}.ITI);
        for k=1:11        
            DLL1(subj,k) = Data{subj}.SimSerie{1}.fit{k}.LL; %/ntaps; % - Data{subj}.SimSerie{1}.fit{1}.LL;
            DLL2(subj,k) = Data{subj}.SimSerie{2}.fit{k}.LL; %/ntaps; % - Data{subj}.SimSerie{2}.fit{1}.LL;
        end
    end
    
    f32 = figure(32);
    clf;
    delta = 0:5:50;   
    hold on
    
    m1 = sum(DLL1);
	m2 = sum(DLL2);

    l(2) = plot(delta,m1,'k-.','linewidth',2);
    l(1) = plot(delta,m2,'r-','linewidth',2);
    
    
    setaxis2(fsize)
	set(gca,'Xtick',[0 20 40]);
    legend(l,{'$b\neq 1$','$b=1$'},'box','off','interpreter','latex','Location','NorthWest')
    
    
    xlabel('refractory time $\Delta$ [ms]','interpreter','latex')
    ylabel('log-likelihood $L$','interpreter','latex')
    xlim([0 50])
    
    print(f32,'-dpdf',[root filenameFig3b '.pdf'])
    savefig(f32,[root filenameFig3b '.fig'])
    
end


% ITI distribution for one subject
if Fig3c == 1    
   f41 = setfigure(33);
   clf;

   theta = Data{sub}.SimSerie{SimSerieHard}.fit{fitHard}.theta;
   PlotITIContDataFit(Data,sub,theta,ModelHard,1,fsize)
   ylim([10^(-15) 1])
   set(gca,'YTick',[1E-15,1E-10,1E-5,1E0])
   
   print(f41,'-dpdf',[root filenameFig3c '.pdf'])
   savefig(f41,[root filenameFig3c '.fig'])
end

% ITI distribution for all subjects
if Fig3d == 1
    f42 = setfigure(34);
    clf;    
    
    for k=1:length(Data)        
	PlotITI(Data{k}.ITI,1,'color',[0.7 0.7 0.7])
        if k==1
            hold on
        end
    end
   
    x = logspace(1,8,100)';
    
    % computes the average fitted ITI
	THETAHard = ExtractField(Data,['SimSerie{' num2str(SimSerieHard) '}.fit{' num2str(fitHard) ' }.theta' ]);
    for k=1:length(Data)
        MPCont(:,k) = PCont(x,THETAHard(k,:)',ModelHard);        
    end
    mPCont = mean(MPCont,2);
    
    loglog(x,mPCont,'color',[0 0 0],'linewidth',2);
    setaxis2(fsize)
    xlabel('ITI $\tau$ [ms]','interpreter','latex')
    ylabel('$p(\tau)$','interpreter','latex')
    xlim([5E1 1E8])
    ylim([10^(-15) 1])
    set(gca,'XTick',[1E2 1E4 1E6 1E8])
    set(gca,'YTick',[1E-15,1E-10,1E-5,1E0])
     
    print(f42,'-dpdf',[root filenameFig3d '.pdf'])
    savefig(f42,[root filenameFig3d '.fig'])        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIGURE 4: Fitting results rel. ref  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any([Fig4a Fig4b Fig4c Fig4d Fig4e Fig4f Fig4g Fig4h ])    
    
    THETA = ExtractField(Data,['SimSerie{' num2str(SimSerie) '}.fit{' num2str(fit) ' }.theta' ]);
    la = THETA(:,1);
    lb = THETA(:,2);
    lalph = THETA(:,3);            
    disp(['median(a) = ' num2str(median(la))])
    disp(['median(b) = ' num2str(median(lb))])
    disp(['median(rho) = ' num2str(median(exp(lalph)))])
end    

% ITI distribution for one subject
if Fig4a == 1    
   f41 = setfigure(41);
   clf;

   theta = Data{sub}.SimSerie{SimSerie}.fit{fit}.theta;
   PlotITIContDataFit(Data,sub,theta,Model,1,fsize)
   ylim([10^(-15) 1])
   set(gca,'YTick',[1E-15,1E-10,1E-5,1E0])
   
   print(f41,'-dpdf',[root filenameFig4a '.pdf'])
   savefig(f41,[root filenameFig4a '.fig'])
end

% ITI distribution for all subjects
if Fig4b == 1
    f42 = setfigure(42);
    clf;    
    
    for k=1:length(Data)        
	PlotITI(Data{k}.ITI,1,'color',[0.7 0.7 0.7])
        if k==1
            hold on
        end
    end
   
    x = logspace(1,8,100)';
    
    % computes the average fitted ITI
    for k=1:length(Data)
        MPCont(:,k) = PCont(x,THETA(k,:)',Model);        
    end
    mPCont = mean(MPCont,2);
    
    loglog(x,mPCont,'color',[0 0 0],'linewidth',2);
    setaxis2(fsize)
    xlabel('ITI $\tau$ [ms]','interpreter','latex')
    ylabel('$p(\tau)$','interpreter','latex')
    xlim([5E1 1E8])
    ylim([10^(-15) 1])
    set(gca,'XTick',[1E2 1E4 1E6 1E8])
    set(gca,'YTick',[1E-15,1E-10,1E-5,1E0])
     
    print(f42,'-dpdf',[root filenameFig4b '.pdf'])
    savefig(f42,[root filenameFig4b '.fig'])        
end


% refractory kernel inferred from the subject
if Fig4c == 1    

   f43 = setfigure(43);
   clf;
   
   theta = Data{sub}.SimSerie{SimSerie}.fit{fit}.theta;
   %lambda = Data{sub}.fit{sim}.lambda;
   
   
   PlotKernel(theta,Model,fsize);
   set(gca,'xtick',[10^2, 10^3])
   ylim([0 1.5])
   xlim([20 5000])
    
   print(f43,'-dpdf',[root filenameFig4c '.pdf'])
   savefig(f43,[root filenameFig4c '.fig'])
end

% refractory kernel for all subjects
if Fig4d == 1    

    f44 = setfigure(44);
    clf;
    
    for k=1:length(Data)
        theta = Data{k}.SimSerie{SimSerie}.fit{fit}.theta;
        PlotKernel(theta,Model,fsize,'color',[0.7 0.7 0.7],'linewidth',2);
        if k == 1
            hold on
        end
        [tau,ker(k,:)] = ComputeKernel(theta,Model);
    end
    set(gca,'xtick',[10^2, 10^3])
    ylim([0 1.5])
    xlim([20 5000])
    
   plot(tau,mean(ker),'k','linewidth',2)  % averaged kernel over all subjects
    
   print(f44,'-dpdf',[root filenameFig4d '.pdf'])
   savefig(f44,[root filenameFig4d '.fig'])
end

% priority distribution inferred for the subject
if Fig4e == 1    

   f45 = setfigure(45);
   clf;
   
   %theta = Data{sub}.fit{sim}.theta;  
   theta = Data{sub}.SimSerie{SimSerie}.fit{fit}.theta;  
   
   pd1 = makedist('Beta','a',theta(1),'b',theta(2));
   x = linspace(0,1,100);   
   plot(x,pdf(pd1,x),'k-','linewidth',2)
   setaxis2(fsize)    
   set(gca,'Xtick',[0 0.5 1])   
   xlabel('priority $x$','interpreter','latex')   
   ylabel('$p(x)$','interpreter','latex')
   
      
   print(f45,'-dpdf',[root filenameFig4e '.pdf'])
   savefig(f45,[root filenameFig4e '.fig'])
end


% priority for all subjects
if Fig4f == 1    

    f52 = setfigure(46);
    clf;
                   
    x = linspace(0,1,100);            
    
    for k=1:length(Data)         
        %theta = Data{k}.fit{sim}.theta;               
        theta = Data{k}.SimSerie{SimSerie}.fit{fit}.theta;               
        pd = makedist('Beta','a',theta(1),'b',theta(2));
        plot(x,pdf(pd,x),'color',[0.7 0.7 0.7],'linewidth',2)    
        if k == 1
            hold on
        end
        pdfk(k,:) = pdf(pd,x);
    end
    ylim([0 3])
        
    plot(x,mean(pdfk),'color',[0 0 0],'linewidth',2)                          
    
   setaxis2(fsize)    
   set(gca,'Xtick',[0 0.5 1])
   xlabel('priority $x$','interpreter','latex')
   ylabel('$p(x)$','interpreter','latex')

    
   print(f52,'-dpdf',[root filenameFig4f '.pdf'])
   savefig(f52,[root filenameFig4f '.fig'])
end

% histogram of parameter a
if Fig4g == 1    
    f53 = setfigure(47);
    clf;
    
    histogram(la,'FaceColor',[0.3 0.3 0.3])
    setaxis2(fsize)    
    xlim([0 1.3])
    
    xlabel('$a$','interpreter','latex')
    ylabel('number of subjects','interpreter','latex')
    
    print(f53,'-dpdf',[root filenameFig4g '.pdf'])
    savefig(f53,[root filenameFig4g '.fig'])
end

% histogram of parameter rho
if Fig4h == 1    
    f53 = setfigure(48);
    clf;
    
    histogram(1000*exp(lalph),12,'FaceColor',[0.3 0.3 0.3])
    setaxis2(fsize)    
    %xlim([0 1.3])
    
    xlabel('$\rho$ [Hz]','interpreter','latex')
    ylabel('number of subjects','interpreter','latex')
    
    print(f53,'-dpdf',[root filenameFig4h '.pdf'])
    savefig(f53,[root filenameFig4h '.fig'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     FIGURE 5: model comparison      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Model comparison with BICtot
if Fig5a == 1
           
   nsub = length(Data);
    ngamma = length(Data{1}.SimSerie{5}.fit); % assuming Serie 5 and 6 have same length
    ndelta = length(Data{1}.SimSerie{1}.fit);
    
    for s=1:nsub
        nITI(s) = length(Data{s}.ITI);
        LL(1,s) = Data{s}.SimSerie{1}.fit{1}.LL;
        LL(2,s) = Data{s}.SimSerie{2}.fit{1}.LL;
        LL(3,s) = Data{s}.SimSerie{1}.fit{end}.LL;
        LL(4,s) = Data{s}.SimSerie{2}.fit{end}.LL;        
        
        for k=1:ngamma
            LL5(k,s) = Data{s}.SimSerie{5}.fit{k}.LL;
            LL6(k,s) = Data{s}.SimSerie{6}.fit{k}.LL;
        end
    end
    
    Ntot = sum(nITI);
    Npar(1) =  Data{1}.SimSerie{1}.fit{1}.Npar;
    Npar(2) =  Data{1}.SimSerie{2}.fit{1}.Npar;
    Npar(3) =  Data{1}.SimSerie{1}.fit{end}.Npar + 1; % since we take the optimal \Delta, it counts as an additional parameter
    Npar(4) =  Data{1}.SimSerie{2}.fit{end}.Npar + 1; % idem
  
	% NOTE1: SimSerie 1 and 2 are simulated with fixed delta (hard refractory time). When delta
	% is fixed, the number of parameters is 2, but one take the maximum
	% likelihood solution across delta, then there are 3 parameters.
	% This is why we take here Npar = 3 instead of 2. The maximum
	% liklihood solution is always when delta = 50 ms and therefore, we
	% take ndelta
	% Note2: this slight change in the computation of BIC makes
	% literally no change in the plot

    
    for m=1:4
        BICtot(m) = log(Ntot)*nsub*Npar(m)-2*sum(LL(m,:));
    end
       
    for k=1:ngamma
        Npar5(k) =  Data{1}.SimSerie{5}.fit{k}.Npar;
        Npar6(k) =  Data{1}.SimSerie{6}.fit{k}.Npar;
        BICtot5(k) = log(Ntot)*nsub*Npar5(k)-2*sum(LL5(k,:));
        BICtot6(k) = log(Ntot)*nsub*Npar6(k)-2*sum(LL6(k,:));
    end
        
	[BICmin5,nstar5] = min(BICtot5);
	[BICmin6,nstar6] = min(BICtot6);
    
    disp(['min BIC M5: ' num2str(BICmin5) ', nstar = ' num2str(nstar5)])
    disp(['min BIC M6: ' num2str(BICmin6) ', nstar = ' num2str(nstar6)])
    disp(['Delta BIC: ' num2str(BICmin6-BICmin5)])
    
    % reference BIC for other distributions: log-normal, Weibul, gamma
    if regen == 1
        for s=1:nsub
            x = Data{s}.ITI;
            [lambda,k,LLW(s)] = WeibullFit(x);
            [lambda,k,LLLN(s)] = LogNormalFit(x);
            [k,theta,LLG(s)] = GammaFit(x);
            disp(['s = ' num2str(s)])
        end
        BICtotW = log(Ntot)*nsub*2-2*sum(LLW);
        BICtotLN =log(Ntot)*nsub*2-2*sum(LLLN);
        BICtotG = log(Ntot)*nsub*2-2*sum(LLG);        
    else        

        BICtotG = 167453787.4919;        
        BICtotW = 156085351.442;
        BICtotLN = 148497958.4634;
    end
    
    
    disp(['BIC Weibull: ' num2str(BICtotW)])
    disp(['BIC log-norma: ' num2str(BICtotLN)])
    disp(['BIC Gamma: ' num2str(BICtotG)])
    
    f51 = figure(51);
    clf;

    ymax = 1.71E8;
    % nonlinear transform of the y-axis
    ybar = 1.5E8; a = 0.1;
    fy = @(y) y.*(y<ybar) + (ybar+a*(y-ybar)).*(y>=ybar);
    %fyinv = @(yp) yp.*(yp<ybar) + (yp/a - ((1-a)/a)*ybar).*(yp>=ybar);
    
    % hard refractoriness, Delta = 0
	hold on
    r = rectangle('Position',[0.2 fy(ybar) 40 fy(ymax)],'FaceColor',[0.95 0.95 0.95],'linestyle','none');
    
    l(1) = plot([1 ngamma],fy(BICtot(1))*[1 1],'k-.','linewidth',2);
    l(2) = plot([1 ngamma],fy(BICtot(2))*[1 1],'r-.','linewidth',2);

    % hard refractoriness
    l(3) = plot([1 ngamma],fy(BICtot(3))*[1 1],'k:','linewidth',2);
    l(4) = plot([1 ngamma],fy(BICtot(4))*[1 1],'r:','linewidth',2);

    % relative refractoriness
    l(5) = plot(1:ngamma,fy(BICtot5),'k','linewidth',2);   
    l(6) = plot(1:ngamma,fy(BICtot6),'r','linewidth',2);       
	
    
    
    %xlim([0 40])
    %ylim([1.44E8 1.49E8])
    %ylim([1.44E8 1.67E8])

	plot(nstar5,fy(BICmin5),'k*','MarkerSize',10)
    plot(nstar6,fy(BICmin6),'r*','MarkerSize',10)
    setaxis2(20)
	xlabel('no. of basis functions $n$','interpreter','latex')
    ylabel('$BIC_{\rm pop} (\times 10^8)$','interpreter','latex')
    
%    text(5,fy(1.45E8),'relative refractoriness','interpreter','latex','FontSize',25)
%    text(5,fy(1.465E8),'hard refractoriness','interpreter','latex','FontSize',25)
%    text(5,fy(1.481E8),'hard refractoriness, $\Delta = 0$','interpreter','latex','FontSize',25)       
    
     text(5,fy(1.45E8),'priority, relative ref.','interpreter','latex','FontSize',25)
     text(5,fy(1.465E8),'priority, hard ref.','interpreter','latex','FontSize',25)
     text(5,fy(1.481E8),'priority, no ref.','interpreter','latex','FontSize',25)       


    text(5,fy(1.689E8),'Gamma','interpreter','latex','FontSize',25)       
    text(5,fy(1.578E8),'Weibull','interpreter','latex','FontSize',25)       
    text(5,fy(1.488E8),'Log-Normal','interpreter','latex','FontSize',25)           
    
   
    lega(1) = plot([1 ngamma],fy(BICtotLN)*[1 1],'b--','linewidth',2);
    lega(2) = plot([1 ngamma],fy(BICtotW)*[1 1],'color',[0 0.5 0],'linestyle','--','linewidth',2);    
    lega(3) = plot([1 ngamma],fy(BICtotG)*[1 1],'color',[1 0.5 0],'linestyle','--','linewidth',2);    
       
    xlim([0 40])    
    ylim(fy([1.44E8 ymax]))
    yval = [1.44 1.46 1.48 1.5 1.6 1.7]*1E8;
    set(gca, 'ytick', fy(yval), 'yticklabel', num2cell(yval/1E8));
    
 
    leg = legend(l(5:6),{'$b=1$','$b\neq 1$'},'interpreter','latex');
    set(leg,'Position',[0.7114    0.125    0.1823    0.1650],'box','off');
	print(f51,'-dpdf',[root filenameFig5a '.pdf'])
    savefig(f51,[root filenameFig5a '.fig'])
    
end
% comparison of power-law exponent a
if Fig5b == 1

    Nf = length(Data{1}.SimSerie{5}.fit);
    for k=1:length(Data)
        for f=1:Nf
            lista5(k,f) = Data{k}.SimSerie{5}.fit{f}.theta(1);
            lista6(k,f) = Data{k}.SimSerie{6}.fit{f}.theta(1);
        end        
    end
    
    f52 = figure(52);
    clf;
    
    hold on
    

    e(1) = errorbar(1:Nf,median(lista5),median(lista5)-prctile(lista5,25),prctile(lista5,75)-median(lista5),'k','linewidth',2);
    e(2) = errorbar(1:Nf,median(lista6),median(lista6)-prctile(lista6,25),prctile(lista6,75)-median(lista6),'r','linewidth',2);        
    
    %xlim([0 Nf])
	xlim([0 40])   
    ylim([0 1])
    plot([1 Nf],[0.5 0.5],'k-.','linewidth',2)
    setaxis2(fsize)
	set(gca,'Xtick',0:10:40);
      
    leg = legend(e,{'$b=1$','$b\neq 1$'},'interpreter','latex');
    set(leg,'Position',[0.7    0.3    0.18    0.16],'box','off');
    
    xlabel('no. of basis functions $n$','interpreter','latex')
    ylabel('$a$','interpreter','latex')
    
    print(f52,'-dpdf',[root filenameFig5b '.pdf'])
    savefig(f52,[root filenameFig5b '.fig'])
    
end

% comparison of effective time constant tau^*
if Fig5c == 1

    Nf = length(Data{1}.SimSerie{5}.fit);
    for k=1:length(Data)
        for f=1:Nf
            listtau5(k,f) =  Data{k}.SimSerie{5}.fit{f}.taus;
            listtau6(k,f) = Data{k}.SimSerie{6}.fit{f}.taus;
        end        
    end
    
    f53 = figure(53);
    clf;
    
    hold on
    
    e(1) = errorbar(1:Nf,median(listtau5),median(listtau5)-prctile(listtau5,25),prctile(listtau5,75)-median(listtau5),'k','linewidth',2);    
    e(2) = errorbar(1:Nf,median(listtau6),median(listtau6)-prctile(listtau6,25),prctile(listtau6,75)-median(listtau6),'r','linewidth',2);    
    
    %xlim([0 Nf])
	xlim([0 40])   
    ylim([0 1500])
    setaxis2(fsize)
    set(gca,'Xtick',0:10:40);
    
    leg = legend(e,{'$b=1$','$b\neq 1$'},'interpreter','latex');
    set(leg,'Position',[0.7    0.75    0.18    0.16],'box','off');
    
    xlabel('no. of basis functions $n$','interpreter','latex')
    ylabel('$\tau^* [ms]$','interpreter','latex')
        
    print(f53,'-dpdf',[root filenameFig5c '.pdf'])
    savefig(f53,[root filenameFig5c '.fig'])
    
end

if any([Fig6a Fig6b Fig6binset Fig6c Fig6d Fig6e Fig6f])

       
    Nf = length(Data{1}.SimSerie{5}.fit);
    for k=1:length(Data)
        for f=1:Nf
            listtau5(k,f) = Data{k}.SimSerie{5}.fit{f}.taus;
            listtau6(k,f) = Data{k}.SimSerie{6}.fit{f}.taus;

            lista5(k,f) = Data{k}.SimSerie{5}.fit{f}.theta(1);
            lista6(k,f) = Data{k}.SimSerie{6}.fit{f}.theta(1);

        end        
    end
    
    
    for ff = 1:40    
        x5 = listtau5(:,ff);
        y5 = lista5(:,ff);
    
        x6 = listtau6(:,ff);
        y6 = lista6(:,ff);

        %remove outliers: has not effect and actually not required
        taumax = 10000; %max(x);
        ind5 = find(x5>taumax);
        ind6 = find(x6>taumax);
        xp5{ff} = x5;xp5{ff}(ind5)=[];
        yp5{ff} = y5;yp5{ff}(ind5)=[];        
        mdl5{ff} = fitlm(xp5{ff},yp5{ff}); 
        Rsquared5(ff) = mdl5{ff}.Rsquared.Ordinary;
        slope5(ff) = table2array(mdl5{ff}.Coefficients(2,1));
      
        xp6{ff} = x6;xp6{ff}(ind6)=[];
        yp6{ff} = y6;yp6{ff}(ind6)=[];        
        mdl6{ff} = fitlm(xp6{ff},yp6{ff}); 
        Rsquared6(ff) = mdl6{ff}.Rsquared.Ordinary;
        slope6(ff) = table2array(mdl6{ff}.Coefficients(2,1));
    
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIG 6 scale free "a" and tau %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a vs tau
if Fig6a == 1
    f61 = figure(61);
    clf;
%    f61.Resize = 'off';
    hold on
    
    ff5 = 20;
    ff6 = 21;
    hold on
    plot(xp5{ff5},yp5{ff5},'k.','MarkerSize',12)    
    plot(xp6{ff6},yp6{ff6},'r.','MarkerSize',12)    
    
    taumax1 = 1500;
	nx = linspace(1,taumax1,100);    
	d5 = table2array(mdl5{ff5}.Coefficients);
    d6 = table2array(mdl6{ff6}.Coefficients);
    ny5 = d5(1,1)+d5(2,1)*nx;
    ny6 = d6(1,1)+d6(2,1)*nx;
    p(1) = plot(nx,ny5,'k','linewidth',2);
    p(2) = plot(nx,ny6,'r','linewidth',2);   
            
    ylim([0 1])
    xlim([0 taumax1])
    
	setaxis2(fsize)    
        
    xlabel('$\tau^* [ms]$','interpreter','latex')
    ylabel('$a$','interpreter','latex')
    
	leg = legend(p,{'$b=1$','$b\neq 1$'},'interpreter','latex');
    set(leg,'Position',[0.25    0.27    0.18    0.16],'box','off');

      
    print(f61,'-dpdf',[root filenameFig6a '.pdf'])
    savefig(f61,[root filenameFig6a '.fig'])
end


% R2 as a function of n    
if Fig6b == 1    
    f62 = figure(62);
    clf;
%    f62.Resize = 'off';
    hold on
    
    p(1) = plot(1:40,Rsquared5,'k','linewidth',2);
    p(2) = plot(1:40,Rsquared6,'r','linewidth',2);
        
    setaxis2(fsize)    
    set(gca,'Xtick',0:10:40) 
    xlabel('no. of basis functions $n$','interpreter','latex')
    ylabel('$R^2$','interpreter','latex')
    
	leg = legend(p,{'$b=1$','$b\neq 1$'},'interpreter','latex');
    set(leg,'Position',[0.25    0.27    0.18    0.16],'box','off');
      
    print(f62,'-dpdf',[root filenameFig6b '.pdf'])
    savefig(f62,[root filenameFig6b '.fig'])
end

% R2 inset
if Fig6binset == 1
    
	f621 = figure(621);
    clf;
    %f621.Resize = 'off';
    hold on
    
    x = 2:40;
    c = 1; %1E4;
    plot(x,c*slope5(x),'k','linewidth',4)
    plot(x,c*slope6(x),'r','linewidth',4)
    plot([1 40],[0 0],'k-.','linewidth',4)
        
    setaxis2(60)    
    set(gca,'Ytick',[-5E-4*c 0]) 
    ylim([-7E-4*c 0])
    %xlabel('number of basis functions $n$','interpreter','latex')
    xlabel('$n$','interpreter','latex')
    %ylabel('slope ($\times 10^{-4}$)' ,'interpreter','latex')
    ylabel('slope' ,'interpreter','latex')
      
    print(f621,'-dpdf',[root filenameFig6b 'inset.pdf'])
    savefig(f621,[root filenameFig6b 'inset.fig'])
end

% consistency of a for SimSerie5
if Fig6c == 1    
	f63 = figure(63);    
    clf;
%    f63.Resize = 'off';
    hold on
    

	A = lista5(:,10:30);
    [n1,n2] = size(A);
	Ad = diff(A,1,2);
        
    for f=1:n2
        r = randperm(n1,n1)';
        B(:,f) = A(r,f);
    end
    Bd = diff(B,1);
    
    ad = Ad(1:end);
    x = linspace(-0.5,0.5,50);
    %h(1) = histogram(ad,x,'FaceColor',[1 0 0],'Normalization','pdf');    
    bd = Bd(1:end);
    h(2) = histogram(bd,x,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'Normalization','probability');
    h(1) = histogram(ad,x,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'Normalization','probability');
    
    setaxis2(fsize)    

    leg = legend(h,{'fitted','random'},'interpreter','latex');
    set(leg,'Position',[0.65    0.75    0.18    0.16],'box','off');
    
    xlabel('$a_n-a_{n-1}$','interpreter','latex')
    ylabel('$p(a_n-a_{n-1})$','interpreter','latex')
    
	print(f63,'-dpdf',[root filenameFig6c '.pdf'])
    savefig(f63,[root filenameFig6c '.fig'])

    disp(['a: SimSerie 5: std_interfit = ' num2str(std(ad)) ', std_interindividual = ' num2str(std(bd)) ', ratio = ' num2str(std(bd)/std(ad))])
end

% consistency of a for SimSerie6
if Fig6d == 1        
    
	f64 = figure(64);    
    clf;
%    f64.Resize = 'off';
    hold on
    

	A = lista6(:,10:30);
    [n1,n2] = size(A);
	Ad = diff(A,1,2);
        
    for f=1:n2
        r = randperm(n1,n1)';
        B(:,f) = A(r,f);
    end
    Bd = diff(B,1);
    
    ad = Ad(1:end);
    x = linspace(-0.5,0.5,50);   
    bd = Bd(1:end);
    h(2) = histogram(bd,x,'FaceColor',[1 1 1],'EdgeColor',[1 0 0],'Normalization','probability');
    h(1) = histogram(ad,x,'FaceColor',[1 0 0],'EdgeColor',[1 0 0],'Normalization','probability');
    
    setaxis2(fsize)    
    
    leg = legend(h,{'fitted','random'},'interpreter','latex');
    set(leg,'Position',[0.65    0.75    0.18    0.16],'box','off');
    
    xlabel('$a_n-a_{n-1}$','interpreter','latex')
    ylabel('$p(a_n-a_{n-1})$','interpreter','latex')
    
	print(f64,'-dpdf',[root filenameFig6d '.pdf'])
    savefig(f64,[root filenameFig6d '.fig'])
    
    disp(['a: SimSerie 6: std_interfit = ' num2str(std(ad)) ', std_interindividual = ' num2str(std(bd)) ', ratio = ' num2str(std(bd)/std(ad))])
end

% SimSerie5 consistency in taus
if Fig6e == 1
    
    f65 = figure(65);
    clf;
%    f65.Resize = 'off';
    hold on
    

	A = listtau5(:,10:30);
    [n1,n2] = size(A);
	Ad = diff(A,1,2);
        
    for f=1:n2
        r = randperm(n1,n1)';
        B(:,f) = A(r,f);
    end
    Bd = diff(B,1);
    
    ad = Ad(1:end);
    x = linspace(-1000,1000,50);    
    bd = Bd(1:end);
    h(2) = histogram(bd,x,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'Normalization','probability');
    h(1) = histogram(ad,x,'FaceColor',[0 0 0],'EdgeColor',[0 0 0],'Normalization','probability');
    
    setaxis2(fsize)    
    set(gca,'Xtick',-1000:500:1000)
    leg = legend(h,{'fitted','random'},'interpreter','latex');
    set(leg,'Position',[0.65    0.75    0.18    0.16],'box','off');
    
    xlabel('$\tau^*_n-\tau^*_{n-1}$  [ms]','interpreter','latex')
    ylabel('$p(\tau^*_n-\tau^*_{n-1})$','interpreter','latex')
    
	print(f65,'-dpdf',[root filenameFig6e '.pdf'])
    savefig(f65,[root filenameFig6e '.fig'])

    disp(['tau: SimSerie 5: std_interfit = ' num2str(std(ad)) ', std_interindividual = ' num2str(std(bd)) ', ratio = ' num2str(std(bd)/std(ad))])
end

        
% SimSerie6 consistency in taus    
if Fig6f == 1
    f66 = figure(66);
    clf;
    %f66.Resize = 'off';
    hold on
    

	A = listtau6(:,10:30);
    [n1,n2] = size(A);
	Ad = diff(A,1,2);
        
    for f=1:n2
        r = randperm(n1,n1)';
        B(:,f) = A(r,f);
    end
    Bd = diff(B,1);
    
    ad = Ad(1:end);
    x = linspace(-1000,1000,50);
    bd = Bd(1:end);
    h(2) = histogram(bd,x,'FaceColor',[1 1 1],'EdgeColor',[1 0 0],'Normalization','probability');
    h(1) = histogram(ad,x,'FaceColor',[1 0 0],'EdgeColor',[1 0 0],'Normalization','probability');
    
    setaxis2(fsize)    
    set(gca,'Xtick',-1000:500:1000)
	leg = legend(h,{'fitted','random'},'interpreter','latex');
    set(leg,'Position',[0.65    0.75    0.18    0.16],'box','off');
    
    xlabel('$\tau^*_n-\tau^*_{n-1}$ [ms]','interpreter','latex')
    ylabel('$p(\tau^*_n-\tau^*_{n-1})$','interpreter','latex')
    
	print(f66,'-dpdf',[root filenameFig6f '.pdf'])
    savefig(f66,[root filenameFig6f '.fig'])

    disp(['tau: SimSerie 6: std_interfit = ' num2str(std(ad)) ', std_interindividual = ' num2str(std(bd)) ', ratio = ' num2str(std(bd)/std(ad))])        
    
    
end





