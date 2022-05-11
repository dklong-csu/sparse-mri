%% This script will calculate the SVD for each saved forward model

clear;clc;
% Get the name of every saved model
fileInfo = dir('./matrix-data/*.mat');
fnames = {fileInfo.name};
dirs = {fileInfo.folder};

% iterate through each file, load the matrix, perform SVD, and export the image
for i=1:length(fnames)
    load(fullfile(dirs{i},fnames{i}),"dft2D_mtx_sparse_real")
    S = svd(dft2D_mtx_sparse_real);
    figure
    plot(S)
    split_file = split(fnames{i},'.');
    fileName = strcat('svd_',split_file{1},'.png');
    filePath = fullfile('pics',fileName);
    exportgraphics(gca, filePath,'Resolution',800)
end