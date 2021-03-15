%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

function [rangeDataPeak,xM,yM] = phaseProfileArbitraryArray2D(sarDataIn,iParams,fParams)
%% Declare Optional Parameters
%-------------------------------------------------------------------------%
if iParams.CSAR
    iParams.nUsefulHorMeasurement = iParams.nAngMeasurement;
end

%% Resize the Data
%-------------------------------------------------------------------------%
if ismatrix(sarDataIn) && min(size(sarDataIn) == [iParams.nVerMeasurement*iParams.nVx,fParams.adcSample])
    % MIMO-SAR 1D
    % sarDataIn = s(y*nVx,k)
    sarDataIn = permute(sarDataIn,[2,1]);
    % sarDataIn = s(k,y*nVx)
else
    warning("Input SAR data is not properly sized")
    return;
end

if ~ismatrix(sarDataIn)
    warning("Input SAR data does not have enough dimensions")
end

sarDataFFT = fft(sarDataIn,iParams.nFFT,1);

%% Find Recommended Peak for Phase
%-------------------------------------------------------------------------%
tempRange = mag2db(squeeze(abs(sarDataFFT(:,:))));
maxIdx = zeros(1,size(tempRange,2));

for ii = 1:size(tempRange,2)
    [~,maxIdx(ii)] = max(tempRange(:,ii));
end

maxIdx = mode(maxIdx);

rangeAxis = generateRangeAxis(fParams,iParams.nFFT);

%% Get User Input for Range Peak to Use
%-------------------------------------------------------------------------%
rangeBinIdx = input("What is the range bin index? (Recommended: " + maxIdx + " -> " + rangeAxis(maxIdx) +  "m) ");

if isempty(rangeBinIdx)
    rangeBinIdx = 0;
end

try
    rangeDataPeak = squeeze(sarDataFFT(rangeBinIdx,:));
catch
    warning("Invalid input, using recommended")
    rangeDataPeak = squeeze(sarDataFFT(maxIdx,:,:));
end

phaseData = unwrap(angle(rangeDataPeak));

%% Display the Phase Profile
%-------------------------------------------------------------------------%
yM = (0:iParams.nVerMeasurement*iParams.nVx-1)*iParams.yStepM_mm;

figure
plot(yM,phaseData)
title("Phase Profile at Range Bin: " + maxIdx + " -> " + rangeAxis(maxIdx) +  "m")

%% Find Phase Minimum Point
%-------------------------------------------------------------------------%
[~,idxVer] = min(phaseData(:));

disp("Minimum point at y index: " + idxVer + " or y = " + yM(idxVer))

end