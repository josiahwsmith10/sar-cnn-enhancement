%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

function [sarDataMIMO,sarDataSISO] = SAR_1D_createEcho_MIMO_SISO(iParams,fParams,p)
%% Declare p(y,k): sarData
%-------------------------------------------------------------------------%
sarDataMIMO = complex(single(zeros(iParams.nVx*iParams.nVerMeasurement,fParams.adcSample)));
sarDataSISO = sarDataMIMO;

%% Declare Wavenumber Vector
%-------------------------------------------------------------------------%
f0 = fParams.f0 + fParams.ADCStartTime*fParams.K; % This is for ADC sampling offset
f = single(f0 + (0:fParams.adcSample-1)*fParams.K/fParams.fS); % wideband frequency

c = 299792458; % physconst('lightspeed'); in m/s
k = 2*pi*f/c;
k = reshape(k,1,[]);
clear f f0

%% Create Virtual Synthetic Aperture Axis
%-------------------------------------------------------------------------%
yV_m = single(zeros(1,iParams.nVx*iParams.nVerMeasurement));

for indVer = 1:iParams.nVerMeasurement
    yV_m( ((indVer-1)*iParams.nVx + 1):(indVer*iParams.nVx) ) = iParams.V_m(:,2) + (indVer-1)*iParams.yStepTot_mm*1e-3;
end

yVmean_m = mean(yV_m);
yV_m = reshape(yV_m - yVmean_m,[],1);

%% Define p(y,z): pyz
%-------------------------------------------------------------------------%
if mod(p.xLim,2) == 0
    pyz = single(squeeze(p.pxyz(end/2,:,:)));
else
    pyz = single(squeeze(p.pxyz((end+1)/2,:,:)));
end

% but if the parallel pool is not open
for iyP = 1:p.yLim
    for izP = 1:p.zLim
        if pyz(iyP,izP) > 1e-8
            for indVer = 1:iParams.nVerMeasurement
                for indTx = 1:iParams.nTx
                    yR_m = iParams.locRx_m(:,2) + (indVer-1)*iParams.yStepTot_mm*1e-3 - yVmean_m;
                    yT_m = iParams.locTx_m(indTx,2)  + (indVer-1)*iParams.yStepTot_mm*1e-3 - yVmean_m;
                    Rr = sqrt( (p.yT(iyP) - yR_m).^2 + (p.zT(izP)).^2 );
                    Rt = sqrt( (p.yT(iyP) - yT_m).^2 + (p.zT(izP)).^2 );
                    idxY = ((indVer-1)*iParams.nVx + (indTx-1)*iParams.nRx + 1):((indVer-1)*iParams.nVx + indTx*iParams.nRx);
                    sarDataMIMO(idxY,:) = sarDataMIMO(idxY,:) + pyz(iyP,izP) .* (Rt.*Rr).^(-1) .* exp(1j*k.*(Rt + Rr) - 1j*pi*fParams.K/(c^2).*(Rr + Rr).^2);
                end
            end
            R = sqrt( (p.yT(iyP) - yV_m).^2 + (p.zT(izP)).^2 );
            sarDataSISO = sarDataSISO + R.^(-2) .* exp(1j*2*k.*R);
        end
    end
end

%% Display the Reflectivity Function
%-------------------------------------------------------------------------%
figure;
mesh(p.zT,p.yT,pyz,'FaceColor','interp','LineStyle','none');
title("Reflectivity Function of Target Scene")
xlabel("z (m)")
ylabel("y (m)")