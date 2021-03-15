clear;

%% Add Paths
%-------------------------------------------------------------------------%
addpath(genpath("../../../../sar-simulation-jws"))

%% Create p
%-------------------------------------------------------------------------%
clear p
p.xLim = 256;
p.zLim = 256;

p.xMax = 0.15;
p.zMax = 0.15;

p.xT = linspace(-p.xMax+2*p.xMax/p.xLim,p.xMax,p.xLim);
p.zT = linspace(-p.zMax+2*p.zMax/p.zLim,p.zMax,p.zLim);

p.pxz = zeros(p.xLim,p.zLim);

%% Load fParams and iParams
%-------------------------------------------------------------------------%

load fParamsAll; load iParamsAll; load pAll
fParams = fParamsAll.v3;                    % Frequency Parameters
iParams = iParamsAll.SISO_CSAR;             % Image and Scanning Parameters (360deg of rotation)
clear fParamsAll iParamsAll pAll

iParams.nAngMeasurement = 1000;
iParams.tStepM_deg = 360/iParams.nAngMeasurement;
iParams.R0_mm = 250;
iParams.nFFT = 2048;
iParams.xU = 1;
iParams.zU = 1;
iParams.displayResult = false;
iParams.PFA = "linear";

%% Directory for Files
%-------------------------------------------------------------------------% 
ds_path = "D:/Research/SAR Image Despeckling/CSAR 1D/Huge/";

%% Generate Training Data GPU (MOST EFFICIENT)
%-------------------------------------------------------------------------%
numSample = 50000;
numTargetMax = 200;

iParams.displayRandResult = false;

tic
f0 = fParams.f0 + fParams.ADCStartTime*fParams.K; % This is for ADC sampling offset
f = f0 + (0:fParams.adcSample-1)*fParams.K/fParams.fS; % wideband frequency

c = 299792458; % physconst('lightspeed'); in m/s
k = 2*pi*f/c;
k = gpuArray(single(reshape(k,1,[])));
clear f f0
theta_rad = ( (-iParams.nAngMeasurement/2):( (iParams.nAngMeasurement/2) - 1) )*iParams.tStepM_deg*pi/180;
theta_rad = gpuArray(single(reshape(theta_rad,[],1)));

[X,Z] = ndgrid(p.xT,p.zT);

% kX and kZ
kX = 2*k.*cos(theta_rad);
kZ = 2*k.*sin(theta_rad);

kXmax = max(max(kX));
kZmax = max(max(kZ));
kXmin = min(min(kX));
kZmin = min(min(kZ));

kXU = reshape(linspace(iParams.xU*kXmin,iParams.xU*kXmax,iParams.xU*iParams.nFFT),[],1);
kZU = reshape(linspace(iParams.zU*kZmin,iParams.zU*kZmax,iParams.zU*iParams.nFFT),1,[]);

kXUsteps = mean(diff(kXU));
kZUsteps = mean(diff(kZU));

% Upsample kU and theta_radU
theta_radU = atan2(kZU,kXU);
kU = (1/2)*sqrt(kXU.^2 + kZU.^2);
kUpsample = reshape(linspace(min(k),max(k),length(kZU)),1,[]);
theta_radUpsample = reshape(linspace(min(theta_rad),max(theta_rad),length(kXU)),[],1);

% Azimuth Filter
azimuthFilter = exp(1j*2*kUpsample*(iParams.R0_mm*1e-3).*cos(theta_radUpsample));
azimuthFilterFFT = fft(azimuthFilter,[],1);

% Region of Interest
xRangeT_m = (-(iParams.nFFT*iParams.xU)/2:((iParams.nFFT*iParams.xU)/2-1))*(2*pi/(kXUsteps*length(kXU)));
zRangeT_m = (-(iParams.nFFT*iParams.zU)/2:((iParams.nFFT*iParams.zU)/2-1))*(2*pi/(kZUsteps*length(kZU)));
indX = xRangeT_m >= (p.xT(1)) & xRangeT_m <= (p.xT(end));
indZ = zRangeT_m >= (p.zT(1)) & zRangeT_m <= (p.zT(end));

% Input and Output
Input = gpuArray(single(zeros(p.xLim,p.zLim)));
Output = gpuArray(single(zeros(p.xLim,p.zLim)));

for indSample = 1:numSample
    % Input and Output
    csarData = gpuArray(single(complex(zeros(iParams.nAngMeasurement, fParams.adcSample))));
    
    % Target Scene
    numTarget = randi(numTargetMax,1);
    coord_x_m = p.xT(1) + (p.xT(end)-p.xT(1))*rand(numTarget,1);
    coord_z_m = p.zT(1) + (p.zT(end)-p.zT(1))*rand(numTarget,1);
    amp_xz = 0.5 + 0.5*rand(numTarget,1);
    if indSample == 1
        numTarget = 1;
        coord_x_m = p.xT(end/2);
        coord_z_m = p.zT(end/2);
        amp_xz = 1;
    end
    R0 = iParams.R0_mm*1e-3;
    for indTarget = 1:numTarget
        R = sqrt( (R0*cos(theta_rad) - coord_x_m(indTarget)).^2 + (R0*sin(theta_rad) - coord_z_m(indTarget)).^2 );
        csarData = csarData + amp_xz(indTarget) .* (R.^(-2)) .* exp(1j*2*k.*R - 1j*pi*fParams.K*(2*R/c).^2);
        Output = Output + abs(amp_xz(indTarget))*exp( - (X - coord_x_m(indTarget)).^2/0.001^2 - (Z - coord_z_m(indTarget)).^2/0.001^2);
    end
    
    csarDataUpsampled = interpn(theta_rad,k,csarData,theta_radUpsample,kUpsample,iParams.PFA,0);
    csarDataUpsampledFFT = fft(csarDataUpsampled,[],1);
    azimuthFiltered = conj(ifft(csarDataUpsampledFFT .* conj(azimuthFilterFFT),[],1));
    
    csarImageFFT = interpn(theta_radUpsample,kUpsample,azimuthFiltered,theta_radU,kU,iParams.PFA,0);
    
    csarImage = abs(ifftshift(ifft2(csarImageFFT)));
    
    csarImage = csarImage(indX,indZ);
    
    Input = imresize(csarImage,[p.xLim,p.zLim]);
    
    % Normalize Input and Output
    Input = gather(Input./max(abs(Input(:))));
    Output = gather(Output./max(abs(Output(:))));
    
    % Save the Input and Output Images
    save(ds_path + "Input/Input_" + indSample + ".mat","Input")
    save(ds_path + "Output/Output_" + indSample + ".mat","Output")
end
clear indSample amp_xz azimuthFilter azimuthFilterFFT azimuthFiltered c coord_x_m coord_z_m csarDataUpsampled csarDataUpsampledFFT csarData csarImage csarImageFFT indTarget indX indZ k kU kUpsample kX kXmax kXmin kXU kXUsteps kZ kZmax kZmin kZU kZUsteps R R0 theta_rad theta_radU theta_radUpsample X xRangeT_m Z zRangeT_m
toc 

%% Save the Training Data

save("./Training Data/trainingSet","Input","Output","numTargetMax","numSample")