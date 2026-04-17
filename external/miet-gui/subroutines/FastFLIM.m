function FastFLIM(tag,tim,thresh)

figure

ind = tag<thresh;
tim(ind) = 0;
tim = max(tim,0);

ind = ~ind(:);
tmp = sort(tim(ind));
tmp = round(tmp*10)./10;
tmp = tmp(tmp>0);
n = round(numel(tmp)/50);
lims = [tmp(n) tmp(49*n)];

intens = tag.^0.8;
intens = intens./max(intens(:));

xl = [399 440 485 540 580 610 630 699];
yl = [[0 0 0]; [0 0 1]; [0 1 1]; [0 1 0]; [1 1 0]; [1 0.65 0]; [1 0 0]; [1 0 1]];
lambda   = 399:3:699;
spectrum = interp1(xl, yl, lambda);
colormap(spectrum)

val = tim(:);

val(isnan(val))  = 0;
val(val<lims(1)) = 0;
val(val>lims(2)) = lims(2);

k      = 2 + round((val-lims(1))./(lims(2)-lims(1)).*(size(spectrum,1)-2));
k(k<2) = 1;

im = spectrum(k,:);
im = reshape(im,size(tag,1),size(tag,2),3);
im = im.*repmat(intens, [1 1 3]);

image(im)

set(gca,'DataAspectRatio', [1,1,1], ...
    'PlotBoxAspectRatio',[1 1 1], ...
    'XDir','normal', ...
    'YDir','reverse', ...
    'FontSize',9,...
    'CLim', lims, ...
    'color', [0 0 0]);

n= 1;
for ttst = lims(1):0.5:lims(2)
    stst(n) = cellstr(sprintf('%4.1f',ttst));
    n = n+1;
end
colorbar('Limits', lims, 'YTickMode','manual','XColor',[0 0 0],'YColor',[0 0 0],'Box','off','YTick', lims(1):0.5:lims(2), ...
    'FontSize', 9,'TickDirection','out', 'YTickLabel',stst);

