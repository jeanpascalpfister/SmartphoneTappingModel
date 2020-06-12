function varargout = PlotITI(x,varargin)
% JPP 23.9.2015
% 31.5.2017 added varargout

% plots the ITI intervals in loglog
if nargin>1
    scale = varargin{1};
else
    scale = 1;
end

if nargin>2
    opt = varargin(2:end);
else
    opt{1} = 'ko';
end

% if nargin>2
%     opt = varargin{:};
% else
%     opt = 'ko';
% end

f= @(x) p*min(max(0.5+c*x,0),1);


xmind = min(x(x~=-Inf));
xmax = max(x(x~=Inf));

bins = logspace(log10(xmind),log10(xmax),50);
[N,edges] = histcounts(x,bins,'Normalization','pdf');
e = edges(1:end-1);        %TBC

NN = N;
%NN = N/sum(N);

%figure;
%clf;
h = loglog(scale*e,NN,opt{:});
setaxis2(25)
%xlim([1 10^8])
%ylim([10^(-15) 1])

if nargout > 0
    varargout{1} = h;
end
