function [f1,f2] = PlotDBIC(BIC,modelref,model,threshold,varargin)
% JPP 13.3.2018
% Plots the histogram of the difference in BIC 

if nargin>4
    n1 = varargin{1};
    f1 = setfigure(n1);
else
    f1 = setfigure;
end


colL = [1 0 0]; % color for DBIC < threshold
colM = [0 0 0]; % color for DBIC < threshold       
colH = [0 0 1]; % color for DBIC > -threshold

db = BIC(:,modelref)-BIC(:,model); % \Delta BIC

% computes the bin edges such that one edge is at zero
r = abs(min(db))/(max(db)-min(db));
NbL = round(r*20);
xL = linspace(min(db),0,NbL);
dx = xL(2)-xL(1);
x = [xL,dx:dx:(max(db)+dx)];

histogram(db(db<threshold),x,'FaceColor',colL);
hold on    
%histogram(db(abs(db)<abs(threshold)),x,'FaceColor',colM);
histogram(db(db>-threshold),x,'FaceColor',colH);
a = axis;
plot([0 0],a(3:4),'k-.','linewidth',2);

setaxis2(25);
ylabel('number of subjects','interpreter','latex')
str = ['$BIC(M_' num2str(modelref) ')-BIC(M_' num2str(model) ')$'];
xlabel(str,'interpreter','latex')


% inset
if nargin>5
    n2 = varargin{2};
    f2 = setfigure(n2);
else
    f2 = setfigure;
end

xmax = 100;
Nbins = 21;
x = linspace(-xmax,xmax,Nbins);
dbinset = db; dbinset = dbinset(abs(dbinset)<xmax);
histogram(dbinset(dbinset<threshold),x,'FaceColor',colL)
hold on
histogram(dbinset(dbinset>-threshold),x,'FaceColor',colH)
histogram(dbinset(abs(dbinset)<abs(threshold)),x,'FaceColor',colM)
setaxis2(40);
plot(threshold*ones(1,2),ylim,'k-.')
plot(-threshold*ones(1,2),ylim,'k-.')
ylabel('n. of subjects','interpreter','latex')
xlabel(str,'interpreter','latex')
xlim([-xmax xmax])    
y = ylim;
ylim([0 y(2)]);

end

