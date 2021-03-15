%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

%% Clear
%-------------------------------------------------------------------------%
clear
addpath(genpath("./"))

%% Create Datastore
%-------------------------------------------------------------------------%

inputData = fileDatastore(fullfile(ds_path + "Input/"),"ReadFcn",@load,"FileExtensions",".mat");
outputData = fileDatastore(fullfile(ds_path + "Output/"),"ReadFcn",@load,"FileExtensions",".mat");

inputData = transform(inputData,@(data) inputPreprocess(data));
outputData = transform(outputData,@(data) outputPreprocess(data));
ds = combine(inputData,outputData);

%% Save Datastore
%-------------------------------------------------------------------------%
save smallMIMO ds
