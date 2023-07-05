%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

%% 1. Add the Necessary Folders to Path (Run First)
%-------------------------------------------------------------------------%
addpath(genpath("../../"),genpath("../../../sar-simulation-jws/"))

%% Josiah's Custom Huge MIMO
%-------------------------------------------------------------------------%
lambda_m = 299792458/(79e9);

% Rx Antenna 1 is the reference
% Coordinates: [x y] x-Horizontal, y-Vertical

iParams.locRx_m =  [0       0
                    0       1
                    0       2
                    0       3
                    0       4
                    0       5
                    0       6
                    0       7
                    0       8
                    0       9
                    0       10
                    0       11
                    0       12
                    0       13
                    0       14
                    0       15]*lambda_m/2;

iParams.TxRxOffsetx_m = 0;
iParams.TxRxOffsety_m = 15*lambda_m/2+5e-3;

iParams.locTx_m =  [0       0
                    0       8
                    0       16
                    0       24
                    0       32
                    0       40
                    0       48
                    0       56]*lambda_m + [iParams.TxRxOffsetx_m,iParams.TxRxOffsety_m];

% Choose Active Antennas and Order

iParams = arbitraryUniformArrayDimensions(iParams);

%% Create iParams
%-------------------------------------------------------------------------%
iParams.nVerMeasurement = 4;
iParams.z0_mm = 500;
iParams.nFFT = iParams.nVerMeasurement*iParams.nVx;
iParams.lambda_mm = physconst('lightspeed')/(79e9)*1e3;
iParams.scanName = "Big MIMO"

%% Create 3D Grid of Point Reflectors
%-------------------------------------------------------------------------%
clear p

p.xLim = 256;
p.yLim = 256;
p.zLim = 256;

p.xMax = 0.15;
p.yMax = 0.25;
p.zMax = 0.5;

p.xT = linspace(-p.xMax+2*p.xMax/p.xLim,p.xMax,p.xLim);
p.yT = linspace(-p.yMax+2*p.yMax/p.yLim,p.yMax,p.yLim);
p.zT = linspace(2*p.zMax/p.zLim,2*p.zMax,p.zLim);

p.pxyz = zeros([length(p.xT),length(p.yT),length(p.zT)]);
p.pxyz((end/2-16):16:(end/2+16),(end/2-16):16:(end/2+16),(end/2-32):32:(end/2+32)) = 1;

%% fParams: v3
%-------------------------------------------------------------------------%
fParams.K = 19.988e12;
fParams.fS = 322e3;
fParams.adcSample = 64;
fParams.TXStartTime = 1e-6;
fParams.ADCStartTime = 6e-6;
fParams.RampEndTime = 200.12e-6;
fParams.f0 = 77e9;

%% 3. Get the MIMO Echo Signal sarData = s(x,k)
%-------------------------------------------------------------------------%
% Scanning is performed in the horizontal domain only using MIMO
% transceiver
[sarDataYMIMO,sarDataYSISO] = SAR_1D_createEcho_MIMO_SISO(iParams,fParams,p);

%% 4. Perform Phase Correction
%-------------------------------------------------------------------------%
sarDataYPC = phaseCorrectionArbitraryArray(sarDataYMIMO,iParams,fParams);

%% 5. Reconstruct the 2D Image using the Range Migration Algorithm (RMA)
%-------------------------------------------------------------------------%
sarImageY_2D_RMA_MIMO = SAR_1D_reconstructImage_2D_RMA_MIMO(sarDataYPC,iParams,fParams);

%% 5. Reconstruct the 2D Image using the Range Migration Algorithm (RMA)
%-------------------------------------------------------------------------%
sarImageY_2D_RMA_SISO = SAR_1D_reconstructImage_2D_RMA_MIMO(sarDataYSISO,iParams,fParams);

%% 6. Reconstruct the 2D Image using the Back Projection Algorithm (BPA)

sarImageY_2D_BPA = SAR_1D_reconstructImage_2D_BPA_MIMO(sarDataYPC,iParams,fParams,p);

%% 7. Reconstruct the 2D Image using the Matched Filter (MF) Technique

sarImageY_2D_MF = SAR_1D_reconstructImage_2D_MF_MIMO(sarDataYPC,iParams,fParams,p);