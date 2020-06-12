function [f1,f2,f3,f4] = PlotNiceTaps(varargin)
% JPP 3.11.2016
% plot taps at various zooms 

% import data
if nargin > 0
    Data = varargin{1};
else
    fileFitExp = '../results/DataFitExp.mat';
    load(fileFitExp,'DataFitExp');    
    Data = FilterData(DataFitExp); 
end

%%%%%%%%%%%%%%%%%%%%%%%
% parameters of fig 1 %
%%%%%%%%%%%%%%%%%%%%%%%

k = 1; % subject number
v1 = [2015,04,18];
v2 = [2015,04,19];
unit = 1/(3600*1000); % number of ms per unit
unitstr = 'hr';
xt = [0 6 12 18 24];

box = [0, 800, 800, 250];

f1 = setfigure(1,box);
clf;
D = Data{k};
taps = (ExtractTaps(D,v1,v2) - TimeVec2Stamp(v1))*unit;
PlotTaps(taps,unitstr)
xlim([min(xt) max(xt)])
set(gca,'xtick',xt)

[daynum,dayname] = weekday(datenum(v1),'long');
disp(dayname);
%title([dayname ' ' datestr(datenum(v1),29)],'Fontsize',20)


%%%%%%%%%%%%%%%%%%%%%%%
% parameters of fig 2 %
%%%%%%%%%%%%%%%%%%%%%%%

k = 1; % subject number
v1 = [2015,04,18,6,0,0];
v2 = [2015,04,18,7,0,0];
unit = 1/(60*1000); % number of ms per unit
unitstr = 'min';
xt = 0:10:60;


f2 = setfigure(2,box);
clf;
D = Data{k};
taps = (ExtractTaps(D,v1,v2) - TimeVec2Stamp(v1))*unit;
PlotTaps(taps,unitstr)
xlim([min(xt) max(xt)])
set(gca,'xtick',xt)

%%%%%%%%%%%%%%%%%%%%%%%
% parameters of fig 3 %
%%%%%%%%%%%%%%%%%%%%%%%

k = 1; % subject number
v1 = [2015,04,18,6,2,0];
v2 = [2015,04,18,6,3,0];
unit = 1/(1000); % number of ms per unit
unitstr = 's';
xt = 0:10:60;


f3 = setfigure(3,box);
clf;
D = Data{k};
taps = (ExtractTaps(D,v1,v2) - TimeVec2Stamp(v1))*unit;
PlotTaps(taps,unitstr)
xlim([min(xt) max(xt)])
set(gca,'xtick',xt)


%%%%%%%%%%%%%%%%%%%%%%%
% parameters of fig 4 %
%%%%%%%%%%%%%%%%%%%%%%%

k = 1; % subject number
v1 = [2015,04,18,6,2,45];
v2 = [2015,04,18,6,2,51];
unit = 1/(1000); % number of ms per unit
unitstr = 's';
xt = 45:51;


f4 = setfigure(4,box);
clf;
D = Data{k};
taps = 45 + (ExtractTaps(D,v1,v2) - TimeVec2Stamp(v1))*unit;
PlotTaps(taps,unitstr)
xlim([min(xt) max(xt)])
set(gca,'xtick',xt)


end

