function PlotTaps(taps,varargin)
% JPP 3.11.2016
%
% PlotTaps(taps(1:100))

if nargin>1
    unitstr = varargin{1};
else
    unitstr = 'ms';
end

hold on;

n = length(taps);
for k=1:n
    t = taps(k);
    plot([t,t],[0,1],'k','linewidth',2)
end

setaxis2(50);
xlabel(['time [' unitstr ']'],'interpreter','latex')
set(gca,'YTickLabel',[])
set(gca,'YTick',[])
set(gca,'Ycolor',[1 1 1])
ylim([0 1.5])

end

