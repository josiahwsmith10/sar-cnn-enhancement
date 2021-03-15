%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

function [Ry_m,Ty_m,Vy_m,Dy_m] = linearMIMOArrayDimensions(nTx,delyTx,nRx,delyRx,TxRxOffset_m,plotArrays)
%
% Rx (o) Spacing - lambda/2
% Tx (x) Spacing - 2*lambda
%
%                               -> positive y direction (for scanner)
%          Rx      Tx       Tx
%       o o o o     x        x


%% Define Optional Parameters
%-------------------------------------------------------------------------%
if nargin < 6
    plotArrays = false;
end

%% Declare Spatial Features of Physical Array
%-------------------------------------------------------------------------%
Ry_m = (0:(nRx-1)) * delyRx;                                     % Location of physical receive elements
Ty_m = TxRxOffset_m + (nRx-1) * delyRx + (0:(nTx-1)) * delyTx;    % Location of physical transmit elements

%% Declare Difference Vector
%-------------------------------------------------------------------------%
Dy_m = (Ty_m' - Ry_m).';
Dy_m = Dy_m(:).'; % Distance between each Tx/Rx pair

%% Declare Spatial Features of Virtual Array
%-------------------------------------------------------------------------%
Vy_m = (Ry_m + Ty_m.').'/2;
Vy_m = Vy_m(:).'; % Location of virtual elements

%% Plot the Arrays
%-------------------------------------------------------------------------%
if plotArrays
    figure
    scatter(Ry_m,zeros(length(Ry_m),1))
    hold on
    scatter(Ty_m,zeros(length(Ty_m),1),'x')
    
    scatter(Vy_m,ones(size(Vy_m)))
    ylim([-0.5 1.5])
    xlabel("y (m)")
    legend("Rx","Tx","Virtual Array")
end

%% Check for Overlapping Elements
%-------------------------------------------------------------------------%
if size(Vy_m,2) ~= size(unique(Vy_m),2)
    warning(size(Vy_m,2) - size(unique(Vy_m),2) + " elements are overlapping, out of " + size(Vy_m,2))
end