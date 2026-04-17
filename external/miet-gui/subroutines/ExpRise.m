function [err, c, zz, z] = ExpRise(p, t, y, pic, pos)

t = t(:);%-min(abs(t(:)));
p = 1./p(:)';
y = y(:);
% tmp = y(1)./y;
tmp = y;

zz = [ones(numel(t),1) -exp(-abs(t)*p)];

if nargin>5 && ~isempty(pos)
    c = lsqnonneg(zz,tmp);
else
    c = zz\(tmp);
end
z = zz*c;
if nargin>3 && ~isempty(pic)
    if pic==1
%         subplot(4,1,1:3)
        plot(t, y, 'ob', t, z, 'r');
%         axis tight
        ylabel('excitation chance')
%         subplot(4,1,4)
%         plot(t, (tmp-z)./sqrt(abs(z)),t,0*t)
%         axis tight
        xlabel('delay(ns)')
        
        drawnow
    elseif pic==2
        semilogy(t, y, 'o', t, z); drawnow
    else
        semilogx(t, y, 'o', t, z); drawnow
    end
end

% err = sum(sum((tmp-z).^2./abs(z)));
err = sum((tmp-z).^2);
