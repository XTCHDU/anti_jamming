function [valueTOA,valueTOE,vSigData] = TOA(vSignalInput)
%INPUTS:
% TOA---Time of arrival estimation, including time of end estimation.
%vSignalInput
%   vSignalInput means a vector of radar signals after channelization.
%   data format of vSignalInput is usually complex double.
%OUTPUTS:
%valueTOA
%   TOA value, usually expressed as corresponding index. 
%valueTOE
%   Time of end value, usually expressed as corresponding index. 
%vSigData
%   vSigData means a vector of signal samples in the corresponding channel,
%same as vSignalInput.

numEnergy = 8;
numPoints = length(vSignalInput);
vSignalInputTemp = [vSignalInput,zeros(1,numEnergy)];

vSignalFlag = zeros(1,numPoints);


valueThre = 0.3; %   Threshold needs to be changed. For TOA estimation, choose 0.45 and for PW estimation, choose 0.3  
for kCol=1:numPoints
    sigPower = mean(abs(vSignalInputTemp(kCol:kCol + numEnergy)).^2);
    sigPowerTemp(kCol) = sigPower;
    if sigPower >= valueThre
       vSignalFlag(kCol) = 1;
    end    
end

vSigSmooth = smooth(vSignalFlag);   % values belong to{0,0.2,0.4,0.6,0.8,1}
valueThre2 = 0.5;

vFlagTemp = zeros(1,numPoints);
for kCol = 1:numPoints
    if vSigSmooth(kCol) > valueThre2
        vFlagTemp(kCol) = 1;
    end
end
vFlagTemp2 = diff([0,vFlagTemp,0]);

valueTOATemp = find(vFlagTemp2 == 1);
valueTOETemp = find(vFlagTemp2 == -1) - 1;

if isempty(valueTOATemp)
    valueTOA = 1;
else
    valueTOA = valueTOATemp(1) + numEnergy/2;
end
if isempty(valueTOETemp)
    valueTOE = numPoints;
else
    valueTOE = valueTOETemp(length(valueTOETemp)) + numEnergy/2;
end

if valueTOE>length(vSignalInput)
    valueTOE = length(vSignalInput);
end

vSigData = vSignalInput(valueTOA:valueTOE);

