%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

function [netOut,Layers,Options] = fcnn_trainNetwork2D(inputData,outputData)
%% Get Data Dimensions and Create Input Data Matrix
%-------------------------------------------------------------------------%
if ndims(inputData) == 3
    inputData = permute(inputData,[1 2 4 3]);
    outputData = permute(outputData,[1 2 4 3]);
end

if ~isa(inputData,"single")
    inputData = single(gather(inputData));
    outputData = single(gather(outputData));
end

%% Get the Network Layers and Options
%-------------------------------------------------------------------------%
[Layers,Options] = fcnn_createNetwork2D(size(inputData,[1,2,3]));

%% Train the Network
%-------------------------------------------------------------------------%
disp('Training Network!')
netOut = trainNetwork(inputData,outputData,Layers,Options);
end