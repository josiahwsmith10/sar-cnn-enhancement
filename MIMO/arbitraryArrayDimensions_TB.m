%% Copyright(C) 2020 The University of Texas at Dallas
%  Developed by: Josiah W. Smith
%  Advisor: Prof. Murat Torlak
%  Department of Electrical and Computer Engineering

%  Redistributions and use of source must retain the above copyright notice
%  Redistributions in binary form must reproduce the above copyright notice

%% AWR1443
%-------------------------------------------------------------------------%

lambda_m = 299792458/(79e9);

dTxRx = 5e-3;

% Rx Antenna 1 is the reference
% Coordinates: [x y] x-Horizontal, y-Vertical

iParams.locRx_m =  [0   0
                    0   lambda_m/2
                    0   lambda_m
                    0   3*lambda_m/2];

% 12 Channel
iParams.locTx_m =  [0               3*lambda_m/2+dTxRx
                    -lambda_m/2     3*lambda_m/2+dTxRx+lambda_m
                    0               3*lambda_m/2+dTxRx+2*lambda_m   ];
                
% 8 Channel
% iParams.locTx_m =  [0               3*lambda_m/2+dTxRx
%                     0               3*lambda_m/2+dTxRx+2*lambda_m   ];

iParams = arbitraryUniformArrayDimensions(iParams,true);

%% Cascaded AWR2243
%-------------------------------------------------------------------------%

lambda_m = 299792458/(79e9);

% Rx Antenna 1 is the reference
% Coordinates: [x y] x-Horizontal, y-Vertical

iParams.locRx_m =  [0       0       
                    0       0.5     
                    0       1       
                    0       1.5     
                    0       19.5    
                    0       20      
                    0       20.5    
                    0       21      
                    0       17.5    
                    0       18      
                    0       18.5    
                    0       19      
                    0       -5.5    
                    0       -5      
                    0       -4.5    
                    0       -4]*lambda_m;

iParams.TxRxOffsetx_m = 10e-3;
iParams.TxRxOffsety_m = -20e-3;

iParams.locTx_m =  [0       0      
                    -1      -0.5    
                    -2.5    -1      
                    -3      10.5    
                    -3      8.5     
                    -3      6.5     
                    -3      4.5     
                    -3      2.5     
                    -3      0.5     
                    -3      -1.5    
                    -3      -3.5    
                    -3      -5.5]*lambda_m + [iParams.TxRxOffsetx_m,iParams.TxRxOffsety_m];

% Choose Active Antennas and Order
iParams.locRx_m = iParams.locRx_m([13,14,15,16,1,2,3,4,9,10,11,12,5,6,7,8],:);
iParams.locTx_m = iParams.locTx_m([12,11,10,9,8,7,6,5,4],:);

% Two Device Cascade Mode
% rxAntPos = rxAntPos([4,3,2,1,16,15,14,13],:);
% txAntPos = txAntPos([1,2,3,10,11,12],:);

iParams = arbitraryUniformArrayDimensions(iParams,true);

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

iParams = arbitraryUniformArrayDimensions(iParams,true)