%% This script will perform an iterative reconstruction for each forward model saved

clear;clc;close all;
% Get the name of every saved model
fileInfo = dir('./matrix-data/*.mat');
fnames = {fileInfo.name};
dirs = {fileInfo.folder};

addpath(genpath('./flexBox/'));

reg_parameters = [0.01, 0.1, 1];
% iterate through each file, load the matrix, perform reconstruction, and export the figure
for i=1:length(fnames)
    load(fullfile(dirs{i},fnames{i}))
    split_file = split(fnames{i},'.');
    ssim_file = fullfile("matrix-data",strcat("ssim_",split_file{1},".txt"));
    fileID = fopen(ssim_file,'a+');
    
    for RegParam=reg_parameters
        % Reconstruction and denoising using ROF model
        main = flexBox;
        main.params.tryCPP = 0; %change, if C++ module is compiled

        %add primal var u
        numberU = main.addPrimalVar(size(image));
        
        %add data-fidelity: 1/2\|Au-f\|_2^2
        main.addTerm(L2dataTermOperator(1,dft2D_mtx_sparse_real,b_sparse_real_Noise),numberU);
        
        %add regularizer
        main.addTerm(L1gradientIso(RegParam,size(image)),numberU); %TV-regularization
        iter_type = 'TV';
%         main.addTerm(L2identity(RegParam,size(image)),numberU); % Tikhonov L2 regularization
%         iter_type = 'Tikhonov';
        %main.addTerm(L2laplace(RegParam,size(image)),numberU); % Generalized Tikhonov with Laplace operator
        % iter_type = 'Gen_Tikhonov';
        
        %box constraint ensures the result to stay in [0,1]
        main.addTerm(boxConstraint(0,1,size(image)),numberU);
        
        %run minimization algorithm
        tic;main.runAlgorithm;toc;
        
        %get result
        result = main.getPrimal(numberU);
        
        %show result
        figure
        reconstr_image = reshape(result,[M,M]);
        plotim_final = [image,reconstr_image];
        imagesc(plotim_final)
        axis image
        colormap(gray)
        xticks('')
        yticks('')
        similarity = ssim(reconstr_image, image);

        fprintf(fileID,"%s reconstruction alpha = %f: SSIM=%f\n",iter_type,RegParam,similarity);
        
        fileName = strcat('iterative_',split_file{1},'_',iter_type,'_',num2str(RegParam),'.png');
        filePath = fullfile('pics',fileName);
        exportgraphics(gca, filePath,'Resolution',800)
        close all
        clear main
    end
    fclose(fileID);
end