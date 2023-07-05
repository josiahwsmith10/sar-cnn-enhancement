%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

function iParams = arbitraryUniformArrayDimensions(iParams,plotArrays)
%% Define Optional Parameters
%-------------------------------------------------------------------------%
if nargin < 2
    plotArrays = false;
end

%% Change Locations of Tx and Rx to nTx x nRx x 2
%-------------------------------------------------------------------------%
iParams.nTx = size(iParams.locTx_m,1);
iParams.nRx = size(iParams.locRx_m,1);

txT = reshape(iParams.locTx_m,iParams.nTx,1,[]);
rxT = reshape(iParams.locRx_m,1,iParams.nRx,[]);

%% Declare Difference Vector
%-------------------------------------------------------------------------%
iParams.D_m = reshape(permute(txT-rxT,[2,1,3]),[],2);   % Distance between each Tx/Rx pair

%% Declare Spatial Features of Virtual Array
%-------------------------------------------------------------------------%
iParams.V_m = reshape(permute((txT+rxT)/2,[2,1,3]),[],2);

%% Check for Overlapping Elements
%-------------------------------------------------------------------------%
[~,iParams.idxUnique] = uniquetol(iParams.V_m,'Byrows',true);
if size(iParams.V_m,1) ~= size(uniquetol(iParams.V_m,'Byrows',true),1)
    warning(size(iParams.V_m,1) - size(uniquetol(iParams.V_m,'Byrows',true),1) + " elements are overlapping, out of " + size(iParams.V_m,1))
    iParams.V_m = uniquetol(iParams.V_m,'Byrows',true);
    iParams.D_m = uniquetol(iParams.D_m,'Byrows',true);
end
iParams.nVx = numel(iParams.idxUnique);

%% Set yStepM_mm and xStepM_mm
%-------------------------------------------------------------------------%
if numel(uniquetol(iParams.V_m(:,2))) ~= 1
    iParams.yStepM_mm = mean(diff(iParams.V_m(:,2)))*1e3;
    iParams.yStepTot_mm = size(iParams.V_m,1)*iParams.yStepM_mm;
else
    warning("Not setting iParams.yStepM_mm: all elements are at the same height")
end
if numel(uniquetol(iParams.V_m(:,1))) ~= 1
    iParams.xStepM_mm = mean(diff(iParams.V_m(:,1)))*1e3;
    iParams.xStepTot_mm = size(iParams.V_m,2)*iParams.xStepM_mm;
else
    warning("Not setting iParams.xStepM_mm: all elements are at the same horizontal position")
end

%% Plot the Arrays
%-------------------------------------------------------------------------%
if plotArrays
    figure
    scatter(iParams.locRx_m(:,1)*1e3,iParams.locRx_m(:,2)*1e3)
    k=1:length(iParams.locRx_m(:,1)); text(iParams.locRx_m(:,1)*1e3+0.05,iParams.locRx_m(:,2)*1e3+0.05,num2str(k'));
    hold on
    scatter(iParams.locTx_m(:,1)*1e3,iParams.locTx_m(:,2)*1e3,'x')
    k=1:length(iParams.locTx_m(:,1)); text(iParams.locTx_m(:,1)*1e3+0.05,iParams.locTx_m(:,2)*1e3+0.05,num2str(k'));
    xlabel("x (mm)")
    ylabel("y (mm)")
    legend("Rx","Tx")
    
    figure
    scatter(iParams.V_m(:,1)*1e3,iParams.V_m(:,2)*1e3)
    k=1:length(iParams.V_m(:,1)); text(iParams.V_m(:,1)*1e3+0.05,iParams.V_m(:,2)*1e3+0.05,num2str(k'));
    xlabel("x (mm)")
    ylabel("y (mm)")
    legend("Virtual Array")
    
end