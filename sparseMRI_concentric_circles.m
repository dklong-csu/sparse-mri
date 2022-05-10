%% Preliminaries

% Set RNG for reproducible results
rng(0,'twister')

% Graphical parameters
fsize = 24;
dividercolor = .98; % Grayscale of vertical bars between sub-images. Zero=black, one=white
dividerwidth = .05; % Width of vertical bars between sub-images, relative to phantom width
gammacorr = .7; % Adjusting image brightness. Between 0 and 1 will brighten, over 1 will darken

% Read in the MRI image
im = imread('pics/SamuBrain_256.png');
im = double(im); % Go from uint8 type to floating point type
im = im/255; % Normalize image pixel values between 0 and 1
[row,col] = size(im);

% Apply FFT
Fim = fft2(im);

vert_divider = dividercolor*ones(row,round(dividerwidth*col));

%% Reconstruction
Fim_sparse5 = zeros(size(Fim));

% Create index vector indicating the known spokes, and form index image

t = linspace(-1,1,row); % Assuming row=col
[X5,Y5] = meshgrid(t); % Coordinates for the (fftshift-rearranged) frequency domain

nr = 33; %number of concentric circles

maxr = 0.5; %maximum allowed additional radius away from inner circle

fullr = 0.0; %radius of the cicrle we may choose to always include

indim5 = zeros(size(Fim));

for i = 1:nr
%select the appropriate radius
R = i*maxr/nr + fullr;%<<uncomment this to add the circles outside of the inner one 

indim5(abs(abs(X5+1i*Y5)-R*ones(row))<0.01) = 1; % form the circle

end

% option to always include the center circle.
 indim5(abs(X5 + 1i*Y5)<fullr) = 1;

index5 = (fftshift(indim5)>0);

% Apply the index vector to frequency domain information
Fim_sparse5(index5) = Fim(index5);
% Apply inverse Fourier transform
filtered_im5 = real(ifft2(Fim_sparse5));
% Scale filtered image to interval [0,1]. Note that after this step the
% pixel values are not anymore directly comparable to the original image
% pixel values.
filtered_im5 = max(0,filtered_im5); % Remove possible negative pixels
filtered_im5 = filtered_im5/max(filtered_im5(:)); % Scale max to 1

ssim_index5 = ssim(filtered_im5,im);

compression5 = size(find(indim5==1),1)/(row*col);

% Take a look
figure
clf
plotim5 = [im,vert_divider,indim5,vert_divider,filtered_im5];
imagesc(plotim5.^gammacorr)
axis equal
axis off
colormap gray
title(['Concentric Circles ' , 'ssim: ', num2str(ssim_index5) , ' comp: ' , num2str(compression5)],'fontsize',9)

% Save image to file
% imwrite(uint8(255*plotim5.^gammacorr),'pics/MRI_sparsedemo_5_6.png','png')