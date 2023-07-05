function [csarData,pxz] = CSAR_1D_createEcho_SISO_randPoints(iParams,fParams,p)

%% Define Optional Parameters
%-------------------------------------------------------------------------%
if ~isfield(iParams,"displayRandResult")
    iParams.displayRandResult = false;
end

%% Maintain everything in the size of:
% nAngMeasurement x adcSample
%-------------------------------------------------------------------------%
csarData = single(complex(zeros(iParams.nAngMeasurement, fParams.adcSample)));

%% Declare Wavenumber Vector
%-------------------------------------------------------------------------%
f0 = fParams.f0 + fParams.ADCStartTime*fParams.K; % This is for ADC sampling offset
f = f0 + (0:fParams.adcSample-1)*fParams.K/fParams.fS; % wideband frequency

c = 299792458; % physconst('lightspeed'); in m/s
k = 2*pi*f/c;
k = single(reshape(k,1,[]));
clear f f0

%% Declare theta Synthetic Aperture Vector
%-------------------------------------------------------------------------%
theta_rad = ( (-iParams.nAngMeasurement/2):( (iParams.nAngMeasurement/2) - 1) )*iParams.tStepM_deg*pi/180;
theta_rad = single(reshape(theta_rad,[],1));

%% Create Target Scene
%-------------------------------------------------------------------------%
coord_x_m = p.xT(1) + (p.xT(end)-p.xT(1))*rand(p.numTarget,1);
coord_z_m = p.zT(1) + (p.zT(end)-p.zT(1))*rand(p.numTarget,1);

amp_xz = abs(randn(p.numTarget,1));

pxz = single(zeros(p.xLim,p.zLim));

[X,Z] = ndgrid(p.xT,p.zT);

%% Create Echo Signal
%-------------------------------------------------------------------------%
R0 = iParams.R0_mm*1e-3; % m
if ~isempty(gcp('nocreate')) % if parallel pool is open
    parfor indTarget = 1:p.numTarget
        R = sqrt( (R0*cos(theta_rad) - coord_x_m(indTarget)).^2 + (R0*sin(theta_rad) - coord_z_m(indTarget)).^2 );
        csarData = csarData + amp_xz(indTarget) .* (R.^(-2)) .* exp(1j*2*k.*R - 1j*pi*fParams.K*(2*R/c).^2);
        pxz = pxz + abs(amp_xz(indTarget))*exp( - (X - coord_x_m(indTarget)).^2/0.001^2 - (Z - coord_z_m(indTarget)).^2/0.001^2);
    end
else
    for indTarget = 1:p.numTarget
        R = sqrt( (R0*cos(theta_rad) - coord_x_m(indTarget)).^2 + (R0*sin(theta_rad) - coord_z_m(indTarget)).^2 );
        csarData = csarData + amp_xz(indTarget) .* (R.^(-2)) .* exp(1j*2*k.*R - 1j*pi*fParams.K*(2*R/c).^2);
        pxz = pxz + abs(amp_xz(indTarget))*exp( - (X - coord_x_m(indTarget)).^2/0.001^2 - (Z - coord_z_m(indTarget)).^2/0.001^2);
    end
end

if iParams.displayRandResult
    figure;
    subplot(121)
    mesh(p.zT,p.xT,pxz,'FaceColor','interp','LineStyle','none');
    title("Reflectivity Function of Target Scene")
    xlabel("z (m)")
    ylabel("x (m)")
    xlim([p.zT(1) p.zT(end)])
    ylim([p.xT(1) p.xT(end)])
    view(2)
    iParams.displayResult = false;
    csarImage = CSAR_1D_reconstructImage_2D_PFA(csarData,iParams,fParams,p);
    subplot(122)
    mesh(p.zT,p.xT,csarImage,'FaceColor','interp','LineStyle','none');
    title("Reflectivity Function of Target Scene")
    xlabel("z (m)")
    ylabel("x (m)")
    xlim([p.zT(1) p.zT(end)])
    ylim([p.xT(1) p.xT(end)])
    view(2)
end

end