%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

%% Add Necessary Paths
%-------------------------------------------------------------------------%
addpath(genpath("../../"))

%% Load the Datastore
%-------------------------------------------------------------------------%
load("./Training Data/smallMIMO.mat")

%% Train the FCNN
%-------------------------------------------------------------------------%
CSAR_2D_FCNN = fcnn_trainNetwork2D_DS_Custom(ds);

%% Test the FCNN
%-------------------------------------------------------------------------%
indSample = 1;

Input = ds.preview{1};
GroundT = ds.preview{2};

OutputTest = predict(CSAR_2D_FCNN,Input);

figure; mesh(p.zT,p.yT,GroundT); title("Ground Truth Image")

iParams.scanName = "Input";
sarImageInput = SAR_1D_reconstructImage_2D_RMA_MIMO(complex(Input(:,1:8:end,1),Input(:,1:8:end,2)),iParams,fParams,p); view(-90,0)

figure; mesh(p.zT,p.yT,OutputTest); title("Direct Image!")

%% Functions
%-------------------------------------------------------------------------%
function [netOut,Layers,Options] = fcnn_trainNetwork2D_DS_Custom(ds)
%% Get the Network Layers and Options
%-------------------------------------------------------------------------%
[Layers,Options] = fcnn_createNetwork2D_Custom(size(ds.preview{1},[1,2,3]));

%% Train the Network
%-------------------------------------------------------------------------%
disp('Training Network!')
netOut = trainNetwork(ds,Layers,Options);
end

function [Layers,Options] = fcnn_createNetwork2D_Custom(inputSize)
%% Declare Network Layers
%-------------------------------------------------------------------------%
% Have worked for me before
Layers = [ ...
    imageInputLayer(inputSize)
    % 512 x 512 x 2
    convolution2dLayer(64,16)
    tanhLayer
    % 449 x 449 x 16
    convolution2dLayer(64,16)
    tanhLayer
    % 386 x 386 x 16
    convolution2dLayer(64,16)
    tanhLayer
    % 323 x 323 x 16
    convolution2dLayer(64,16)
    tanhLayer
    % 260 x 260 x 16
    convolution2dLayer(5,1)
    % 256 x 256
    regressionLayer
    ];

%% Declare Network Options
%-------------------------------------------------------------------------%
% Options = trainingOptions('sgdm', ...
%     'Momentum',0.8,...
%     'MaxEpochs',25,...
%     'MiniBatchSize',32,...
%     'InitialLearnRate',1e-3, ...
%     'LearnRateSchedule','piecewise', ...
%     'LearnRateDropFactor',0.5, ...
%     'LearnRateDropPeriod',1, ...
%     'executionEnvironment','gpu', ...
%     'Verbose',1, ...
%     'Shuffle','every-epoch', ...
%     'Plots','training-progress');

Options = trainingOptions("adam", ...
    "MaxEpochs",100, ...
    "InitialLearnRate",1e-3,...
    "MiniBatchSize",4, ...
    "Shuffle","every-epoch", ...
    "Plots","training-progress", ...
    "Verbose",true, ...
    "LearnRateSchedule","piecewise", ...
    "LearnRateDropFactor",1, ...
    "LearnRateDropPeriod",10);

end