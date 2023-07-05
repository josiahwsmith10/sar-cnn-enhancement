clear;

%% Add Paths
%-------------------------------------------------------------------------%
addpath(genpath("../../../../sar-simulation-jws"))

%% Create p
%-------------------------------------------------------------------------%
clear p
p.yLim = 256;
p.zLim = 256;

p.yMax = 0.25;
p.zMax = 0.5;

p.yT = linspace(-p.yMax+2*p.yMax/p.yLim,p.yMax,p.yLim);
p.zT = linspace(2*p.zMax/p.zLim,2*p.zMax,p.zLim);

% fParams: v3
%-------------------------------------------------------------------------%
fParams.K = 19.988e12;
fParams.fS = 322e3;
fParams.adcSample = 64;
fParams.TXStartTime = 1e-6;
fParams.ADCStartTime = 6e-6;
fParams.RampEndTime = 200.12e-6;
fParams.f0 = 77e9;

% Josiah's Custom Huge MIMO
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

% Create iParams
%-------------------------------------------------------------------------%
iParams.nVerMeasurement = 4;
iParams.z0_mm = 500;
iParams.nFFT = iParams.nVerMeasurement*iParams.nVx;
iParams.lambda_mm = physconst('lightspeed')/(79e9)*1e3;
iParams.scanName = "Big MIMO";

%% Directory for Files
%-------------------------------------------------------------------------% 
ds_path = "D:/Research/SAR MIMO Calibration/SAR 2D/Training Data";
ds_path = "./";

%% Generate Training Data GPU (MOST EFFICIENT)
%-------------------------------------------------------------------------%
numSample = 1000;
numTargetMax = 1;

iParams.displayRandResult = false;

tic
% Declare Wavenumber Vector
%-------------------------------------------------------------------------%
f0 = fParams.f0 + fParams.ADCStartTime*fParams.K; % This is for ADC sampling offset
f = single(f0 + (0:fParams.adcSample-1)*fParams.K/fParams.fS); % wideband frequency

c = 299792458; % physconst('lightspeed'); in m/s
k = 2*pi*f/c;
k = reshape(k,1,[]);
clear f f0

% Create Virtual Synthetic Aperture Axis
%-------------------------------------------------------------------------%
yV_m = single(zeros(1,iParams.nVx*iParams.nVerMeasurement));

for indVer = 1:iParams.nVerMeasurement
    yV_m( ((indVer-1)*iParams.nVx + 1):(indVer*iParams.nVx) ) = iParams.V_m(:,2) + (indVer-1)*iParams.yStepTot_mm*1e-3;
end

yVmean_m = mean(yV_m);
yV_m = reshape(yV_m - yVmean_m,[],1);

for indSample = 1:numSample
    numTarget = randi(numTargetMax,1);
    
    % Input and Output
    Input = complex(single(zeros(iParams.nVx*iParams.nVerMeasurement,fParams.adcSample)));
    Output = Input;
    
    % Target Scene
    coord_y_m = p.yT(17) + (p.yT(end-16)-p.yT(17))*rand(numTarget,1);
    coord_z_m = p.zT(17) + (p.zT(end-16)-p.zT(17))*rand(numTarget,1);
    amp_yz = 1.5 + 0.1*rand(numTarget,1);
    
%     [Y,Z] = ndgrid(p.yT,p.zT);
%     Image = reshape(amp_yz,1,1,[]).*exp( - (Y - reshape(coord_y_m,1,1,[])).^2/0.001^2 - (Z - reshape(coord_z_m,1,1,[])).^2/0.001^2);
%     Image = sum(Image,3);
%     max(Image(:))
    
    for indTarget = 1:numTarget
        for indVer = 1:iParams.nVerMeasurement
            for indTx = 1:iParams.nTx
                yR_m = iParams.locRx_m(:,2) + (indVer-1)*iParams.yStepTot_mm*1e-3 - yVmean_m;
                yT_m = iParams.locTx_m(indTx,2)  + (indVer-1)*iParams.yStepTot_mm*1e-3 - yVmean_m;
                Rr = sqrt( (coord_y_m(indTarget) - yR_m).^2 + (coord_z_m(indTarget)).^2 );
                Rt = sqrt( (coord_y_m(indTarget) - yT_m).^2 + (coord_z_m(indTarget)).^2 );
                idxY = ((indVer-1)*iParams.nVx + (indTx-1)*iParams.nRx + 1):((indVer-1)*iParams.nVx + indTx*iParams.nRx);
                Input(idxY,:) = Input(idxY,:) + amp_yz(indTarget) .* (Rt.*Rr).^(-1) .* exp(1j*k.*(Rt + Rr) - 1j*pi*fParams.K/(c^2).*(Rr + Rr).^2);
            end
        end
        R = sqrt( (coord_y_m(indTarget)  - yV_m).^2 + (coord_z_m(indTarget)).^2 );
        Output = Output + R.^(-2) .* exp(1j*2*k.*R);
    end
    
    Input = cat(3,real(Input),imag(Input));
    Output = cat(3,real(Output),imag(Output));
    
    % Save the Input and Output Images
    save(ds_path + "Input/Input_" + indSample + ".mat","Input")
    save(ds_path + "Output/Output_" + indSample + ".mat","Output")
end
toc 

%% Test with Reconstruction Algorithms and Ground Truth
%-------------------------------------------------------------------------%
figure; mesh(p.zT,p.yT,Image,'FaceColor','interp','LineStyle','none'); view(2)
sarImageOutput = SAR_1D_reconstructImage_2D_RMA_MIMO(complex(Output(:,:,1),Output(:,:,2)),iParams,fParams,p);
phaseProfileArbitraryArray2D(complex(Output(:,:,1),Output(:,:,2)),iParams,fParams);

InputPC = phaseCorrectionArbitraryArray(complex(Input(:,:,1),Input(:,:,2)),iParams,fParams);
phaseProfileArbitraryArray2D(InputPC,iParams,fParams);
sarImageInput = SAR_1D_reconstructImage_2D_RMA_MIMO(InputPC,iParams,fParams,p);