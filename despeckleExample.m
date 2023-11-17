function despeckleExample()

% read in image and add speckle to it
imageIn=double(imread('lena_256.tif'));
speckledImage=double(imnoise(uint8(imageIn),'speckle',0.02));

% despeckle the image using novel algorithm
cleanedImage=yu_at4SPEK(101,speckledImage,0.01,80,1,1.2,'SHT','th0+');

% show results
imshow([uint8(speckledImage) uint8(cleanedImage)]);
title('speckled image : cleaned image');
