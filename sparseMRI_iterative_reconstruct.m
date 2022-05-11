%% This script will perform an iterative reconstruction for each forward model saved

clear;clc;close all;
% Get the name of every saved model
% fileInfo = dir('./matrix-data/*.mat');
% fnames = {fileInfo.name};
% dirs = {fileInfo.folder};

% Get specific models
fnames = {"forwardmodel_lowpass5.mat", "forwardmodel_lowpass20.mat"};
dirs = {fullfile('.','matrix-data'), fullfile('.','matrix-data')};

addpath(genpath('./flexBox/'));

reg_parameters = 10.^(-2:1:2);
% iterate through each file, load the matrix, perform reconstruction, and export the figure
for i=1:length(fnames)
    load(fullfile(dirs{i},fnames{i}))
    split_file = split(fnames{i},'.');
    ssim_file = fullfile("matrix-data",strcat("ssim_",split_file{1},".txt"));
    fileID = fopen(ssim_file,'a+');
    
    for RegParam=reg_parameters
        % TV Regularization
        % Reconstruction and denoising using ROF model
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
        resultTV = main.getPrimal(numberU);

        %% Tikhonov regularization
        clear main;
        % Reconstruction and denoising using ROF model
        main = flexBox;
        main.params.tryCPP = 0; %change, if C++ module is compiled

        %add primal var u
        numberU = main.addPrimalVar(size(image));
        
        %add data-fidelity: 1/2\|Au-f\|_2^2
        main.addTerm(L2dataTermOperator(1,dft2D_mtx_sparse_real,b_sparse_real_Noise),numberU);
        
        %add regularizer
        main.addTerm(L2identity(RegParam,size(image)),numberU); % Tikhonov L2 regularization
        
        %box constraint ensures the result to stay in [0,1]
        main.addTerm(boxConstraint(0,1,size(image)),numberU);
        
        %run minimization algorithm
        tic;main.runAlgorithm;toc;
        
        %get result
        resultTik = main.getPrimal(numberU);


        
        %show result
        dividercolor = .98; % Grayscale of vertical bars between sub-images. Zero=black, one=white
        dividerwidth = .05; % Width of vertical bars between sub-images, relative to phantom width
        vert_divider = dividercolor*ones(M,round(dividerwidth*M));

        figure
        reconstr_imageTV = reshape(resultTV,[M,M]);
        reconstr_imageTik = reshape(resultTik,[M,M]);
        plotim_final = [image,vert_divider,plotim, vert_divider,reconstr_imageTV, vert_divider, reconstr_imageTik];
        imagesc(plotim_final)
        axis image
        colormap(gray)
        xticks('')
        yticks('')
        similarityIFFT = ssim(plotim, image);
        similarityTV   = ssim(reconstr_imageTV, image);
        similarityTik  = ssim(reconstr_imageTik, image);

        fprintf(fileID,"-------------------------------\n");
        fprintf(fileID,"Reconstruction using alpha = %f\n",RegParam);
        fprintf(fileID,"-------------------------------\n");
        fprintf(fileID,"IFFT SSIM = %f\n", similarityIFFT);
        fprintf(fileID,"TV SSIM = %f\n", similarityTV);
        fprintf(fileID,"Tik SSIM = %f\n",similarityTik);
        
        fileName = strcat('iterative_',split_file{1},'_',num2str(RegParam),'.png');
        filePath = fullfile('pics',fileName);
        exportgraphics(gca, filePath,'Resolution',800)
        close all
        clear main
    end
    fclose(fileID);
end