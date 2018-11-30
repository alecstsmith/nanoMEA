%Pixel Gradient Analysis 
% 
% 1.two-dimensional convolution 
%   1.1. Gaussian low pass filter 
%   1.2. Sobel horizontal edge-emphasize filter (predefined in Matlab Image Analysis Toolbox) 
%   1.2.1 extract the horizontal edges
%   1.3. Transpose the Sobel filter 
%   1.3.1 extract the vertical edges 
%   1.4. the horizontal and vertical edges are combined 
%   to calculate the gradient magnitude of each pixel in the image
% 2. Thresholding the gradient magnitude to determine the contours of the area of interest  
% 3. Calculate the gradient orientation 
%   3.1. determine the angle of the gradient with respect to the x-axis
%   (Y-axis is oriented downward;The angle range [-?/2, ?/2])
%   3.2. calculate the angle that the orthogonal to the gradient made with
%   x-axis(Y-axis is re-oriented upward;The range [0,180])
%   3.3. The calculated values estimate the principal orientations of
%   individual cells 

% Feel free to modify the threshold value to identify the contours of objets.
% you may use loops to run the analysis on multiple images.
% Written by Hojung Cho, 091007

%find the directory
close all;





% ** Change the directory to where you have your images and this .m file
%cd('C:\Users\Jesse\Dropbox\UW\Lab\Kim Lab\DMD\Confocal\20140602_SarcomereStaining\F-actin');

%call input image


% ** This is the name of your image file
im= imread('crispr_dmd_flat2_60x_6c2.tif');


% ** If the image is red use 1 and so forth
im=double(im(:,:,2)); % 1 --> R, 2--> G, 3 --> B


colormap gray
imagesc(im);

%Using predefined filters, perform 2-D convolution to find the area of
%interest

h=fspecial('gaussian');
im=conv2(im,h,'same');
figure,imagesc(im);%figure 

% use Sobel horizontal filter for y-direction
h=fspecial('sobel');
dimdy=conv2(im,h,'same');

%transpose the matrix for x-direction
h=h';% transpose
dimdx= conv2(im,h,'same');
%figure,imagesc(dimdy),colormap gray;%figure 


%calcualte the gradient magnitude
normsquaregrad= dimdx.^2+dimdy.^2; 

%find the gradient that can be applied to individual cell
Threshold=20000;
masc= (normsquaregrad>Threshold); 
%figure,imagesc(masc.*im);

%angle that the gradient makes with the x-axis, y-axis is oriented downward
%the range is [-pi/2,pi/2]
imangle= masc.*atan(dimdy./dimdx);
%figure,imagesc(imangle);

%angle that the orthogonal to the gradient makes with the x-axis.
%y-axis is oriented upward. range is [0,180]
% This esitmates the principal direction of cells.
imangle=masc.*(pi/2-imangle).*(180/pi);
%figure,imagesc(imangle);
indice=find(masc);% indices where masc= 1
v= imangle(indice);% vector containing the vaules of imangle at locations indices



% ** I use this to shift the hisogram so the max value is at 0deg.
% ** Basically I add or substract a certain angle depending on where the
% max peak appears the first time.
v = v - 0;

shift = v > 90;
shift_v = shift*180;

v = v - shift_v;

M=-85:10:85;
figure, hist(v,M);

aa=hist(v,M);


% ** This is the vector you'll want to save and plot later in Excel or
% SigmaPlot
bb = aa./sum(aa);

figure, plot(M,bb);

%xlswrite('/Users/maca2xdang/Desktop/F-actin Stress Fiber Orientation/aa.xls',aa');