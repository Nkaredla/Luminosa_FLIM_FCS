function [ h, n, mx, sigma, xedges ] = plotHistFit( x, limits, step, unit)
%[ h, n, mx, sigma ] = plotHistFit( x, limits, step, unit)
% plots a 1D histogram of the data in x and a gaussian of mean mx and standard
% deviation sigma, where the latter are calculated in the range of 'limits'
%   step    (optional) stepsize of the bins
%   unit    (optional) the unit of x, used to label the histogram
% 
% OUTPUT
%   h   handle to the histogram figure
%   n   the histogram counts
%   mx  the mean of x in the range of limits
%   
% EXAMPLE
% x = randn(1000,1);
% plotHistFit( x, [-4 4])

x=x(:);
range=max(x)-min(x);
if nargin <2 || isempty(limits)
    limits = [min(x) max(x)];
end
if nargin <4 || isempty(unit)
    unit = '';
end
if nargin <3 || isempty(step)
    h=histogram(x,'BinMethod','fd');
else
    xedges = min(x):step:max(x);
    h=histogram(x,xedges);
end
n = h.Values;
xedges = h.BinEdges;
h=gcf;

mx=mean(x(x>=limits(1) & x<=limits(2)));
sigma=std(x(x>=limits(1) & x<=limits(2)));

xx=linspace(limits(1),limits(2),200);
xlim([floor(min(x)), ceil(max(x))])
hold on
plot(xx,(xedges(2)-xedges(1))*sum(x>=limits(1) & x<=limits(2))*1./sqrt(2*pi*sigma^2)*exp(-(xx-mx).^2./(2*sigma^2)),'k')
hold off
xlabel(['bins [' unit ']'])
ylabel('frequency')
x1 = find(xedges>mx - 2*sigma - 0.2*range,1,'first');
tmp = find(xedges>mx + 2*sigma,1,'first');
if ~isempty(tmp)
    x1(2) =tmp;
else
    x1(2) = length(n);
end
[m] = max(n(min(x1,length(n))));
[~, ind] = min([n(min(x1,length(n)))]);
text(xedges(x1(ind)),max(m,max(n)/2),sprintf('%.1f%s%.1f %s\nN_{total}=%i\nN_{fit}=%i',mx,char(177),sigma,unit,sum(~isnan(x)),sum(x>=limits(1) & x<=limits(2))))


% x1 = round(numel(xcenters)/10);
% x2 = round(mx/range+(numel(xcenters)/10));
% [m] = max(n(x1:x2));
% text(xcenters(x1)+mx,max(m,max(n)/2),sprintf('%.2f%s%.2f %s\nN_{total}=%i\nN_{fit}=%i',mx,char(177),sigma,unit,sum(~isnan(x)),sum(nfit)))

end

