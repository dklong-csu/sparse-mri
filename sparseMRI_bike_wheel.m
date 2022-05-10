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

R = 0.7;

% # spokes

n = 36;

% Initialize the sparse frequency domain information
Fim_sparse6 = zeros(size(Fim));
% Create index vector indicating the known spokes, and form index image
t = linspace(-1,1,row); % Assuming row=col
[X,Y] = meshgrid(t); % Coordinates for the (fftshift-rearranged) frequency domain
% Choose random lines as long as N spokes get chosen
indim6 = zeros(size(Fim)); % Initialize index image
for i = 1:n
    theta = 2*pi/n*i; % direction
    dirind = abs(X*cos(theta)+Y*sin(theta))<.01 & abs(X.^2 + Y.^2 )<R;
    indim6(dirind) = 1;
end

%hub

hubR = 0.25;

indim6(abs(X+1i*Y)<hubR) = 1;

%tire
tireR = 0.8;
indim6(abs(abs(X+1i*Y)-tireR*ones(row))<0.05) = 1;

index6 = (fftshift(indim6)>0);
% Apply the index vector to frequency domain information
Fim_sparse6(index6) = Fim(index6);
% Apply inverse Fourier transform
filtered_im6 = real(ifft2(Fim_sparse6));
% Scale filtered image to interval [0,1]. Note that after this step the
% pixel values are not anymore directly comparable to the original image
% pixel values.
filtered_im6 = max(0,filtered_im6); % Remove possible negative pixels
filtered_im6 = filtered_im6/max(filtered_im6(:)); % Scale max to 1
% Take a look
figure(6)
clf
plotim6 = [im,vert_divider,indim6,vert_divider,filtered_im6];
imagesc(plotim6.^gammacorr)
axis equal
axis off
colormap gray
title('Approach 6: Bike Wheel','fontsize',fsize)

compression6 = size(find(indim6 ==1),1 )/(row*col);

similarity6 = ssim(filtered_im6,im);

fprintf("Bike Wheel ssim: %f compression ratio: %f \n", similarity6, compression6);

% Save image to file
imwrite(uint8(255*plotim6.^gammacorr),'pics/MRI_sparsedemo_6_4.png','png')