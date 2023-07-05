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
load csarHuge

%% Train the FCNN
%-------------------------------------------------------------------------%
CSAR_2D_FCNN = fcnn_trainNetwork2D_despeckle(ds);

%% Test the FCNN
%-------------------------------------------------------------------------%
indSample = 2;

OutputTest = predict(CSAR_2D_FCNN,ds.preview{1});

figure;
subplot(311)
mesh(OutputTest); 
title("Denoised")
subplot(312); 
mesh(ds.preview{1});
title("Input")
subplot(313);
mesh(ds.preview{2});
title("Control")

%% Functions

function [Layers,Options] = fcnn_createNetwork2D_despeckle(inputSize)
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
% Layers = [ ...
%     imageInputLayer(inputSize)
%     
%     convolution2dLayer([3 3],4,"Padding","same")
%     leakyReluLayer
%     
%     convolution2dLayer([5 5],12,"Padding","same")
%     leakyReluLayer
%     
%     convolution2dLayer([3 3],6,"Padding","same")
%     leakyReluLayer
%     
%     convolution2dLayer([1 1],1,"Padding","same")
%     
%     regressionLayer
%     ];


% Same as Gao's paper
Layers = [ ...
    imageInputLayer(inputSize)
    
    convolution2dLayer(30,6,"Padding","same")
    leakyReluLayer
    
    convolution2dLayer(30,12,"Padding","same")
    tanhLayer
    
    convolution2dLayer(30,12,"Padding","same")
    
    convolution2dLayer(30,12,"Padding","same")
        
    convolution2dLayer(1,1,"Padding","same")
    
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
    "MaxEpochs",5, ...
    "InitialLearnRate",1e-2,...
    "MiniBatchSize",32, ...
    "Shuffle","every-epoch", ...
    "Plots","training-progress", ...
    "Verbose",true, ...
    "LearnRateSchedule","piecewise", ...
    "LearnRateDropFactor",0.5, ...
    "LearnRateDropPeriod",10);
end

function [netOut,Layers,Options] = fcnn_trainNetwork2D_despeckle(ds)
%% Get the Network Layers and Options
%-------------------------------------------------------------------------%
[Layers,Options] = fcnn_createNetwork2D_despeckle(size(ds.preview{1},[1,2,3]));

%% Train the Network
%-------------------------------------------------------------------------%
disp('Training Network!')
netOut = trainNetwork(ds,Layers,Options);
end