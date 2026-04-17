% plot intensity and lifetime: 
% in the lifetime image, colour codes lifetime and brightness codes intensity;
% both images are smoothed by a circle of radius 2 pixels
% the numbers in the brackets [] in line 12 are the "colorbar" of life_imm

file='U:\Daja\From Anna\Cells_Atto488_Florian\Images\Image_021_Cell4_PS.mat';

load(file);
figure; mim(tag); colorbar; % plot intensity image

tmpL=life_imm; tmpL(isnan(life_imm))=0; tmpT=tag; tmpT(isnan(life_imm))=0;
figure; mim(mConv2(tmpL.*tmpT,Disk(2))./mConv2(tmpT,Disk(2)),mConv2(tag,Disk(2)),[2.5 4],'v');

tmpL=life2-life1; tmpL(isnan(life1)|isnan(life2))=0; tmpT=tag1; tmpT(isnan(life1)|isnan(life2))=0;
figure; mim(mConv2(tmpL.*tmpT,Disk(1))./mConv2(tmpT,Disk(2)),mConv2(tag1,Disk(1)),[-0.2 0.5],'v');

tmpL=2*life2; tmpL(isnan(life2))=0; tmpT=tag1; tmpT(isnan(life2))=0;
figure; mim(mConv2(tmpL.*tmpT,Disk(1))./mConv2(tmpT,Disk(2)),mConv2(tag1,Disk(1)),[2.5 4],'v');

tmpL=life1; tmpL(isnan(life1))=0; tmpT=tag1; tmpT(isnan(life1))=0;
figure; mim(mConv2(tmpL.*tmpT,Disk(1))./mConv2(tmpT,Disk(1)),mConv2(tag1,Disk(1)),[2.5 4],'v');