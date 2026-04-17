function [img, img_read] = ReadTiff_CP(filename, img_first, img_last)
% ReadTiff(filename, img_first, img_last):
% reads all or parts of a tiff image stack

% if there is no argument, we ask the user to choose a file:
if (nargin == 0 || isempty(filename))
    [filename, pathname] = uigetfile('*.tif;*.stk;*.lsm', 'select image file');
    filename = [ pathname, filename ];
end

if (exist(filename,'file') == 0)
    disp('file not found');
    return;
end
    
if (nargin<=1)  img_first = 1; img_last = inf;   end
if (nargin==2)  img_last = img_first;            end

img = [];

%check for matfile
matfilename = regexprep(filename, '.tif|.stk|.lsm', '.mat');

if (exist(matfilename,'file') > 0)
    % load existing file
    load(matfilename, 'img');
    img_ifds = size(img, 3);
else
    warning off MATLAB:tifflib:TIFFReadDirectory:libraryWarning;
    ti = Tiff(filename, 'r');
    
    % check number of images
    img_ifds = 1;
    while (~ti.lastDirectory())
        ti.nextDirectory();
        img_ifds = img_ifds + 1;
    end
    img_width = ti.getTag('ImageWidth');
    img_length = ti.getTag('ImageLength');
    
    img = zeros(img_length, img_width, img_ifds);
    
    ti.setDirectory(1);
    
    for i=1:img_ifds
        if i > 1
            ti.nextDirectory();
        end
        img(:,:,i) = ti.read();
    end
    ti.close();   
    save(matfilename, 'img');  
end

img_first = min([img_first img_ifds]);
img_last = min([img_last img_ifds]);
img = img(:,:,img_first:img_last);

img_read = size(img,3);
end
