function PlotKernel(theta,Model,fsize,varargin)
% JPP 17.8.2017
% mod 21.8.2018: depends on Model

if nargin>3
    opt = varargin(1:end);
else
    opt{1} = 'k-';
    opt{2} = 'linewidth';
    opt{3} = 2;   
end

[tau,g] = ComputeKernel(theta,Model);

semilogx(tau,g,opt{:});
hold on
semilogx(tau([1,end]),[1 1],'k-.');
setaxis2(fsize);
xlabel('$\tau$ [ms]','interpreter','latex')
ylabel('ref. kernel $r(\tau)$','interpreter','latex')

%if nargout > 0
%    varargout{1} = f;
%end

end

