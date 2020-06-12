function PlotLinITIContDataFit(subject,models)
% JPP 24.8.2017
%
% plots the ITI distribution from a given subject overlaied with a given
% data fit

Data = LoadData;


f = setfigure(1);
clf;
ITI = Data{subject}.ITI;
PlotLinITI(ITI); % Data
hold on
x = linspace(0,2000,100)';% 
%x = logspace(1,8,100)';

%map = colormap('hot');
st = 'LL: ';
for m=1:length(models)
    model = models(m);
    theta = Data{subject}.fit{model}.theta;
    lambda = Data{subject}.fit{model}.lambda;
    k = 1+round((m-1)*40/length(models));
    %h(m) = loglog(x,PCont(x,theta,lambda),'color',map(k,:),'linewidth',2);
    h(m) = plot(x,PCont(x,theta,lambda),'color',[1 0 0],'linewidth',3);
    leg{m} = ['model ' num2str(model)];
    st = [st 'LL' num2str(model) ' = ' num2str(Data{subject}.fit{model}.LL)];
end

disp(st)

%legend(h,leg);
setaxis2(25);
xlabel('ITI [ms]','interpreter','latex')
ylabel('$p(\tau)$','interpreter','latex')
%xlim([50 max(x)])
%set(gca,'Xtick',[1E2 1E4 1E6 1E8])

%PlotKernel(theta,lambda)

% if nargout > 0
%     varargout{1} = f;
% end

end

