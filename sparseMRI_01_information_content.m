% Sparse MRI example using Samu's brain image. The goal is to allow only
% some of the 256x256 FFT elements, or "spokes", to be "measured" 
% (it's all simulated here, actually). Then what is the best sampling 
% strategy for the known spokes? This can be an assignment for students.
%
% Jen and Samu April 2022


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

% Check that the FFT and IFFT work as expected
% figure(1000)
% clf
% imagesc([im,real(ifft2(Fim))])
% axis equal
% axis off
% colormap gray

% Choose how many spokes we can use
ratio = .20; % Ratio of known spokes
N = round(ratio*row*col);


%% Approach #1: random spokes

% Initialize the sparse frequency domain information
Fim_sparse1 = zeros(size(Fim));
% Create index vector indicating the known spokes, and form index image
tmp = randperm(row*col);
index1 = tmp(1:N);
indim1 = zeros(size(Fim));
indim1(index1) = 1;
indim1 = fftshift(indim1); % Rearrange the FFT image so that zero frequency is in the middle
% Apply the index vector to frequency domain information
Fim_sparse1(index1) = Fim(index1);
% Apply inverse Fourier transform
filtered_im1 = real(ifft2(Fim_sparse1));
% Scale filtered image to interval [0,1]. Note that after this step the
% pixel values are not anymore directly comparable to the original image
% pixel values.
filtered_im1 = max(0,filtered_im1); % Remove possible negative pixels
filtered_im1 = filtered_im1/max(filtered_im1(:)); % Scale max to 1
% Take a look
figure
clf
vert_divider = dividercolor*ones(row,round(dividerwidth*col));
plotim1 = [im,vert_divider,indim1,vert_divider,filtered_im1];
imagesc(plotim1.^gammacorr)
axis equal
axis off
colormap gray
title('Approach 1: random sampling','fontsize',fsize)

% Save image to file
imwrite(uint8(255*plotim1.^gammacorr),'pics/MRI_sparsedemo_1.png','png')


%% Approach #2: low-pass

% Initialize the sparse frequency domain information
Fim_sparse2 = zeros(size(Fim));
% Create index vector indicating the known spokes, and form index image
t = linspace(-1,1,row); % Assuming row=col
[X,Y] = meshgrid(t); % Coordinates for the (fftshift-rearranged) frequency domain
R = sqrt(N/(row*col));
indim2 = zeros(size(Fim));
indim2(abs(X+1i*Y)<R) = 1; % Already in the fftshift-rearranged form
index2 = (fftshift(indim2)>0);
% Apply the index vector to frequency domain information
Fim_sparse2(index2) = Fim(index2);
% Apply inverse Fourier transform
filtered_im2 = real(ifft2(Fim_sparse2));
% Scale filtered image to interval [0,1]. Note that after this step the
% pixel values are not anymore directly comparable to the original image
% pixel values.
filtered_im2 = max(0,filtered_im2); % Remove possible negative pixels
filtered_im2 = filtered_im2/max(filtered_im2(:)); % Scale max to 1
% Take a look
figure
clf
plotim2 = [im,vert_divider,indim2,vert_divider,filtered_im2];
imagesc(plotim2.^gammacorr)
axis equal
axis off
colormap gray
title('Approach 2: low-pass filter','fontsize',fsize)

% Save image to file
imwrite(uint8(255*plotim2.^gammacorr),'pics/MRI_sparsedemo_2.png','png')

%% Approach #3: horizontal lines

% Initialize the sparse frequency domain information
Fim_sparse3 = zeros(size(Fim));
% Create index vector indicating the known spokes, and form index image
rowind = 1:round(1/ratio):row; % Pick roughly the right amount of lines
indim3 = zeros(size(Fim));
indim3(rowind,:) = 1;
index3 = (fftshift(indim3)>0);
% Apply the index vector to frequency domain information
Fim_sparse3(index3) = Fim(index3);
% Apply inverse Fourier transform
filtered_im3 = real(ifft2(Fim_sparse3));
% Scale filtered image to interval [0,1]. Note that after this step the
% pixel values are not anymore directly comparable to the original image
% pixel values.
filtered_im3 = max(0,filtered_im3); % Remove possible negative pixels
filtered_im3 = filtered_im3/max(filtered_im3(:)); % Scale max to 1
% Take a look
figure
clf
plotim3 = [im,vert_divider,indim3,vert_divider,filtered_im3];
imagesc(plotim3.^gammacorr)
axis equal
axis off
colormap gray
title('Approach 3: horizontal lines','fontsize',fsize)

% Save image to file
imwrite(uint8(255*plotim3.^gammacorr),'pics/MRI_sparsedemo_3.png','png')


%% Approach #4: random lines through origin

% Initialize the sparse frequency domain information
Fim_sparse4 = zeros(size(Fim));
% Create index vector indicating the known spokes, and form index image
t = linspace(-1,1,row); % Assuming row=col
[X,Y] = meshgrid(t); % Coordinates for the (fftshift-rearranged) frequency domain
% Choose random lines as long as N spokes get chosen
indim4 = zeros(size(Fim)); % Initialize index image
while length(find(indim4(:)>0))<N
    theta = rand*2*pi; % Random direction
    dirind = abs(X*cos(theta)+Y*sin(theta))<.01;
    indim4(dirind) = 1;
end
index4 = (fftshift(indim4)>0);
% Apply the index vector to frequency domain information
Fim_sparse4(index4) = Fim(index4);
% Apply inverse Fourier transform
filtered_im4 = real(ifft2(Fim_sparse4));
% Scale filtered image to interval [0,1]. Note that after this step the
% pixel values are not anymore directly comparable to the original image
% pixel values.
filtered_im4 = max(0,filtered_im4); % Remove possible negative pixels
filtered_im4 = filtered_im4/max(filtered_im4(:)); % Scale max to 1
% Take a look
figure
clf
plotim4 = [im,vert_divider,indim4,vert_divider,filtered_im4];
imagesc(plotim4.^gammacorr)
axis equal
axis off
colormap gray
title('Approach 4: random lines through origin','fontsize',fsize)

% Save image to file
imwrite(uint8(255*plotim4.^gammacorr),'pics/MRI_sparsedemo_4.png','png')


%% Approach #6: random lines through origin, of a fixed radius

%We changed code here.

R = 0.7;

% # spokes

n = 30;

% Initialize the sparse frequency domain information
Fim_sparse6 = zeros(size(Fim));
% Create index vector indicating the known spokes, and form index image
t = linspace(-1,1,row); % Assuming row=col
[X,Y] = meshgrid(t); % Coordinates for the (fftshift-rearranged) frequency domain
% Choose random lines as long as N spokes get chosen
indim6 = zeros(size(Fim)); % Initialize index image
for i = 1:n
    theta = 2*pi/n*i; % Random direction
    dirind = abs(X*cos(theta)+Y*sin(theta))<.01 & abs(X.^2 + Y.^2 )<R;
    indim6(dirind) = 1;
end

%hub

hubR = 0.15;

indim6(abs(X+1i*Y)<hubR) = 1;

%tire
tireR = 0.8;
indim6(abs(abs(X+1i*Y)-tireR*ones(row))<0.03) = 1;

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
imwrite(uint8(255*plotim6.^gammacorr),'pics/MRI_sparsedemo_6.png','png')

