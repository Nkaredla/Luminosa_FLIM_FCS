function imout = shiftlines(imin,direction,points,flagyx)
% imout = shiftlines(imin,direction,points) shifts alternate lines in a
% scan image. 
% imin is the input image
% direction=1 indicates shift every odd line left. 0 -> right (default is 1)
% points is the number of pixels for the shift (default is 1)
% flagyx = 0 if an xy scan is performed. Else the program assumes a yx scan
% and direction = 1 indicates shift every odd vertical line up. By default 
% flagyx = 0

if nargin<4||isempty(flagyx)
    flagyx = 0; % xy scan (shift along the x direction)
else
    imin = imin.'; % yx scan (shift along the y direction)
end
if nargin<3 || isempty(points)
    points = 1;
end
if nargin<2 || isempty(direction)
    direction = 1; % shift left
end
[m,n] = size(imin);

if m/2==round(m/2)
    m = m-1;
    imout = zeros(m,n+points);
else
    imout = zeros(m,n+points);
end

if direction
    imout(1:2:2*floor(m/2)+1,1+points:end)= imin(1:2:2*floor(m/2)+1,:);
    imout(2:2:2*floor(m/2),1:n)= imin(2:2:2*floor(m/2),:);
else
    imout(1:2:2*floor(m/2)+1,1:n)= imin(1:2:2*floor(m/2)+1,:);
    imout(2:2:2*floor(m/2),1+points:end)= imin(2:2:2*floor(m/2),:);
end

imout = imout(:,1+points:end-points);
if flagyx
    imout = imout.';
end
