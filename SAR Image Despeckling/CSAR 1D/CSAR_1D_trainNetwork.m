%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

%% Add Necessary Paths
%-------------------------------------------------------------------------%
addpath(genpath("../../"))

%% Load the Training Data
%-------------------------------------------------------------------------%
load("./Training Data/trainingSet.mat")

%% Train the FCNN
%-------------------------------------------------------------------------%
CSAR_2D_FCNN = fcnn_trainNetwork2D(Input,Output);

%% Test the FCNN
%-------------------------------------------------------------------------%
indSample = 2;

OutputTest = predict(CSAR_2D_FCNN,Input(:,:,indSample));

figure;
subplot(311)
mesh(OutputTest); 
title("Denoised")
subplot(312); 
mesh(Input(:,:,indSample));
title("Input")
subplot(313);
mesh(Output(:,:,indSample));
title("Control")