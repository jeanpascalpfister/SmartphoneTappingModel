function f = setfigure(varargin)
% JPP 9.5.2017
% mod JPP 13.3.2018. no input argument is now also possible

% set nice figure properties
if nargin>0
    n = varargin{1};
    f = figure(n);
else
    f = figure;
end

if nargin>1
    pos = varargin{2};
else
    pos = [840, 500, 600, 450];
end

clf;

set(f,'PaperPositionMode','auto')
set(f, 'Position', pos);      
fig_pos = f.PaperPosition;
f.PaperSize = [fig_pos(3) fig_pos(4)];     

end

