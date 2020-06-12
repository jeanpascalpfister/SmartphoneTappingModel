function PlotITIContDataFit(Data,subject,theta,Model,varargin)
% JPP 24.8.2017
% mod 21.8.2018 only one Model
%
% plots the ITI distribution from a given subject overlaied with a given
% data fit

if nargin>4
    dispfit = varargin{1};
else 
    dispfit = 1;
end

if nargin>5
    fsize = varargin{2};
else 
    fsize = 40;
end


%Data = LoadData;


%f = sefigure(1);
%clf;
ITI = Data{subject}.ITI;
PlotITI(ITI); % Data
hold on
x = logspace(1,8,100)';

if dispfit == 1
    loglog(x,PCont(x,theta,Model),'color','k','linewidth',2);
end

%map = colormap('hot');
%st = 'LL: ';
% for m=1:length(models)
%     model = models(m);
%     theta = Data{subject}.fit{model}.theta;
%     lambda = Data{subject}.fit{model}.lambda;    
%     if dispfit == 1
%         k = 1+round((m-1)*40/length(models));
%         h(m) = loglog(x,PCont(x,theta,lambda),'color',map(k,:),'linewidth',2);
%     end
%     leg{m} = ['model ' num2str(model)];
%     st = [st 'LL' num2str(model) ' = ' num2str(Data{subject}.fit{model}.LL)];
% end

%disp(st)

%legend(h,leg);
setaxis2(fsize);
xlabel('ITI $\tau$ [ms]','interpreter','latex')
ylabel('$p(\tau)$','interpreter','latex')
xlim([50 max(x)])
set(gca,'Xtick',[1E2 1E4 1E6 1E8])

%PlotKernel(theta,lambda)

% if nargout > 0
%     varargout{1} = f;
% end

end

