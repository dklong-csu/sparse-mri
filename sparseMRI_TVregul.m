% Apply Total Variation (TV) regularization to the inverse problem of
% recovering an image from a sparse collection of Fourier transform
% coefficients. The routine sparseMRI_02_forwardmodel_comp.m must be
% computed before this.
%
% TV regularization algoroithm is based on the flexBox toolbox available in
% https://github.com/HendrikMuenster/flexBox
% This is the BibTeX reference: 
% @Article{dirks2015flexbox,
%   Title         = {A Flexible Primal-Dual Toolbox},
%   Author        = {Dirks, Hendrik},
%   Journal       = {ArXiv e-prints},
%   Year          = {2016},
%   Month         = mar,
%   Keywords      = {Mathematics - Optimization and Control, Computer Science - Computer Vision and Pattern Recognition, Computer Science - Mathematical Software, I.4, G.1.6, G.4},
%   Primaryclass  = {math.OC}
% }
%
% Samuli Siltanen and Jennifer Mueller, April 2022



%% Preliminaries

clear all;close all;clc;

% Regularization paramater choice
RegParam = .1;

% Load precomputed data
load data/forwardmodel M N b_sparse_real_Noise dft2D_mtx_sparse_real image

% Graphical parameters
fsize = 20;

<<<<<<< Updated upstream
=======
% Image size will be MxM. M should preferably be a power of 2. 
% It seems that the algorithm has difficulties for M>32.
M = 64;

%read MRI example and downsample
image = im2double(imread('pics/SamuBrain_1024.png'));
image = imresize(image,[M,M]);

%show clean input image
figure(1)
clf
imagesc(image)
axis image
colormap(gray)
title('Input image');

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
title('Sanity check (complex): do we get the original back?');


%% Sparsify the measurement

% Number of known spokes
N = round(spoke_percent*M^2); 

% These construct index1 for each of the three sparsifying methods
% Sparsifying method 1: Pick out vertical lines in the frequency domain
% Nlines = ceil(N/M); % How many lines can we have; the last one possibly only partially
% linestep_max = floor(M/Nlines);
% index1 = zeros(Nlines*M,1); % Initial index vector, may be too long
% for iii = 1:Nlines
%     index1((iii-1)*M+[1:M]) = (iii-1)*linestep_max*M+[1:M];
% end 
% index1 = index1(1:N); % Crop to have the correct amount of spokes in the index1 vector

% Sparsifying method 2: random sampling
tmp = randperm(M^2);  % generate the points for the random spoke locations
index1 = tmp(1:N);

% Sparsifying method 3: low-pass filter
% t = linspace(-1,1,M); 
% [X,Y] = meshgrid(t); % Coordinates for the (fftshift-rearranged into intuitive form) frequency domain
% R = sqrt(N/M^2); % Initial candidate for radius
% index1 = find(abs(X+1i*Y)<R); % This will result in only *approximately* the correct number N of known spokes
% while length(index1)<N
%     R = 1.01*R;
%     index1 = find(abs(X+1i*Y)<R);
% end
% index1 = index1(1:N); % Crop to have the correct amount of spokes in the index1 vector

% Build the system matrix using the index1 vector constructed above
dft2D_mtx_sparse = zeros(N,M^2);
for iii = 1:N
        dft2D_mtx_sparse(iii,:) = dft2D_mtx(index1(iii),:);
end

% show inverse transform, but we need to put the zeros back in
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
title('Inverse transformed sparse image, complex version');


%% Let's build a real-valued model by separating the real and imaginary parts. 

% This trick only works for real-valued images. But that is our case here.
% Separate the real and imaginary parts of the FFT matrix
dft2D_mtx_sparse_real = [real(dft2D_mtx_sparse);imag(dft2D_mtx_sparse)];

% Create data vector and add noise
b_sparse_real = dft2D_mtx_sparse_real*image(:);
b_sparse_real_Noise = b_sparse_real + randn(size(b_sparse_real))*0.01;


% Show inverse transform
% imageTransformNoise2 = b_sparse_real_Noise(1:end/2) + 1i*b_sparse_real_Noise((end/2+1):end);
% plotim = real(ifft2(fftshift(reshape(imageTransformNoise2,[M,M]))));
% plotim = max(0,plotim);
% plotim = plotim/max(plotim(:));
% figure(4)
% clf
% imagesc(plotim.^.5)
% axis image
% colormap(gray)
% title('Inverse transformed sparse image, real version');
>>>>>>> Stashed changes



%% Reconstruction and denoising using ROF model

addpath(genpath('../../flexBox/'));
main = flexBox;
main.params.tryCPP = 0; %change, if C++ module is compiled

%add primal var u
numberU = main.addPrimalVar(size(image));

%add data-fidelity: 1/2\|Au-f\|_2^2
main.addTerm(L2dataTermOperator(1,dft2D_mtx_sparse_real,b_sparse_real_Noise),numberU);

%add regularizer
main.addTerm(L1gradientIso(RegParam,size(image)),numberU); %TV-regularization

%box constraint ensures the result to stay in [0,1]
main.addTerm(boxConstraint(0,1,size(image)),numberU);

%run minimization algorithm
tic;main.runAlgorithm;toc;

%get result
result = main.getPrimal(numberU);

%show result
figure(5)
plotim_final = [image,reshape(result,[M,M])];
imagesc(plotim_final)
axis image
colormap(gray)
title('TV regularized reconstruction from sparse MRI data','fontsize',fsize);
