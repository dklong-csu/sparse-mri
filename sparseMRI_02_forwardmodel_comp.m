% Build a matrix model for the forward problem of sparse MRI. The model is
% not optimised for speed or memory usage; rather it relies on explicitly
% constructed matrices. So it is for demonstration purposes only at low 
% resolutions. 
%
% Samuli Siltanen and Jennifer Mueller, April 2022

%% Preliminaries

clear all;close all;clc;

% Set RNG for reproducible results
rng(0,'twister')

% Regularization paramater choice
RegParam = .1;

% Degree of sparsity in the measurement
spoke_percent  = .05; % eg, 0.2 means 20 percent

% Graphical parameters
fsize = 20;

%% Read in data

% Image size will be MxM. M should preferably be a power of 2. 
% It seems that the algorithm has difficulties for M>32.
M = 32;

%read MRI example and downsample
image = im2double(imread('pics/SamuBrain_1024.png'));
image = imresize(image,[M,M]);


%% Generate matrix representative of Fourier transform using brute force

% First build a complex-valued matrix
dft2D_mtx = zeros(M^2);
for jjj = 1:M^2
   tmpim = zeros(M,M);
   tmpim(jjj) = 1;
   Ftmp = fftshift(fft2(tmpim));  % This maps the low-frequency info that matlab puts in the corners of the square domain to a disk in the center, 
   % which is the intuitive human way to view the fft
   dft2D_mtx(:,jjj) = Ftmp(:);
   if mod(jjj,10*M)==0
       disp([num2str(round(100*jjj/M.^2)),'%'])
   end
end

% Check that the matrix dft2D_mtx really implements 2D Fourier transform
sanityCheck_vec_freq_domain = dft2D_mtx*image(:);
sanityCheck = reshape(sanityCheck_vec_freq_domain,[M,M]);
sanityCheck = real(ifft2(fftshift(sanityCheck)));  % Need to fftshift again so that it toggles back to matlab's way of having the LF info in the corners so that it computes correctly
figure(2)
clf
imagesc([image,sanityCheck])
axis image
colormap(gray)
title('Sanity check: original (left), inverse transformed (right)','fontsize',fsize);


%% Sparsify the measurement

% Number of known spokes
N = round(spoke_percent*M^2); 

% Sparsifying method 1: Pick out vertical lines in the frequency domain
% Nlines = ceil(N/M); % How many lines can we have; the last one possibly only partially
% linestep_max = floor(M/Nlines);
% index1 = zeros(Nlines*M,1); % Initial index vector, may be too long
% for iii = 1:Nlines
%     index1((iii-1)*M+1:M) = (iii-1)*linestep_max*M+1:M;
% end 
% index1 = index1(1:N); % Crop to have the correct amount of spokes in the index1 vector

% Sparsifying method 2: random sampling
% tmp = randperm(M^2);  % generate the points for the random spoke locations
% index1 = tmp(1:N);

% Sparsifying method 3: low-pass filter
t = linspace(-1,1,M); 
[X,Y] = meshgrid(t); % Coordinates for the (fftshift-rearranged into intuitive form) frequency domain
R = sqrt(N/M^2); % Initial candidate for radius
index1 = find(abs(X+1i*Y)<R); % This will result in only *approximately* the correct number N of known spokes
while length(index1)<N
    R = 1.01*R;
    index1 = find(abs(X+1i*Y)<R);
end
index1 = index1(1:N); % Crop to have the correct amount of spokes in the index1 vector

% Build the system matrix using the "index1" vector constructed above
dft2D_mtx_sparse = zeros(N,M^2);
for iii = 1:N
        dft2D_mtx_sparse(iii,:) = dft2D_mtx(index1(iii),:);
end

% Show inverse transform. For that we need to put the zeros back in.
imageTransform_Full = zeros(M,M);
imageTransform_Full(index1) = sanityCheck_vec_freq_domain(index1);
plotim = real(ifft2(fftshift(reshape(imageTransform_Full,[M,M]))));
plotim = max(0,plotim);
plotim = plotim/max(plotim(:));
figure(3)
clf
imagesc(plotim.^.5)
axis image
colormap(gray)
title('Inverse FFT image, unknown spokes replaced by 0','fontsize',fsize);


%% Let's build a real-valued model by separating the real and imaginary parts. 

% This trick only works for real-valued images. But that is our case here!
% Separate the real and imaginary parts of the FFT matrix
dft2D_mtx_sparse_real = [real(dft2D_mtx_sparse);imag(dft2D_mtx_sparse)];

% Create data vector and add noise
b_sparse_real = dft2D_mtx_sparse_real*image(:);
b_sparse_real_Noise = b_sparse_real + randn(size(b_sparse_real))*0.01;

% Save the model to disc
save data/forwardmodel M N b_sparse_real_Noise dft2D_mtx_sparse_real image
