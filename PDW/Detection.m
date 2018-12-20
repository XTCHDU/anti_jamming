function [numFlag] = Detection(vSignalInput)
% Detection---singal detection .
%INPUTS:
%vSignalInput
%signal vector
%OUTPUTS:
%numFlag
%detection flag

numLen = length(vSignalInput);
vSignalFlag = zeros(1,numLen);

numEnergy = 10; %points for energy cumulation
numEffective = 10; %effective points
valueThre = 0.1; %Threshold

for kCol=1:numLen-numEnergy
    sigPower = mean(abs(vSignalInput(kCol:kCol+numEnergy)).^2);
    sigPowerTemp(kCol) = sigPower;
    if sigPower >= valueThre
        vSignalFlag(kCol) = 1;
    end
end

vFlag = sum(vSignalFlag);
if vFlag >= numEffective
    numFlag = 1;
else
    numFlag = 0;
end



