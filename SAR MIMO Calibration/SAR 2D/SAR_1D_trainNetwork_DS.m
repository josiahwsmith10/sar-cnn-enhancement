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
load 64MIMO

%% Train the FCNN
%-------------------------------------------------------------------------%
CSAR_2D_FCNN = fcnn_trainNetwork2D_DS_Custom(ds);

%% Test the FCNN
%-------------------------------------------------------------------------%
indSample = 1;

Input = ds.preview{1};
GoundT = ds.preview{2};

OutputTest = predict(CSAR_2D_FCNN,Input);

iParams.scanName = "Ground Truth";
sarImageGroundT = SAR_1D_reconstructImage_2D_RMA_MIMO(complex(GoundT(:,:,1),GoundT(:,:,2)),iParams,fParams,p); view(-90,0)
iParams.scanName = "Input";
sarImageInput = SAR_1D_reconstructImage_2D_RMA_MIMO(complex(Input(:,:,1),Input(:,:,2)),iParams,fParams,p); view(-90,0)
iParams.scanName = "MIMO Corrected";
sarImageOutput = SAR_1D_reconstructImage_2D_RMA_MIMO(complex(OutputTest(:,:,1),OutputTest(:,:,2)),iParams,fParams,p); view(-90,0)

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
% Layers = [ ...
%     imageInputLayer(inputSize)
%     
%     repmat(convolution2dLayer([5 5],5,"Padding","same"),5,1)    
%     
%     convolution2dLayer([1 1],1,"Padding","same")
%     
%     regressionLayer
%     ];


% Have worked for me before
Layers = [ ...
    imageInputLayer(inputSize)
    
    convolution2dLayer([256 30],16,"Padding","same")
    leakyReluLayer
    
    convolution2dLayer([128 30],16,"Padding","same")
    tanhLayer
    
    convolution2dLayer([64 30],16,"Padding","same")
    
    convolution2dLayer([1 1],2,"Padding","same")
    
    regressionLayer
    ];

% Layers = [ ...
%     imageInputLayer(inputSize)
%     
%     convolution2dLayer([33 5],16,"Padding","same")
% %     leakyReluLayer
%     
%     convolution2dLayer([17 5],16,"Padding","same")
% %     leakyReluLayer
%     
%     convolution2dLayer([5 3],16,"Padding","same")
%     
%     convolution2dLayer([1 1],2,"Padding","same")
%     
%     regressionLayer
%     ];

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
    "MiniBatchSize",32, ...
    "Shuffle","every-epoch", ...
    "Plots","training-progress", ...
    "Verbose",true, ...
    "LearnRateSchedule","piecewise", ...
    "LearnRateDropFactor",1, ...
    "LearnRateDropPeriod",10);

end