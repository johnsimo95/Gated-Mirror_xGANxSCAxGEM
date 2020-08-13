function [stepR,interpTime, interpVoltAll] = ImportLTSpiceSweep(name, timeRes, interpMethod)
%IMPORTLTSPICESWEEP Summary of this function goes here
%   Detailed explanation goes here

fileID = fopen(name,'rt');
formatSpecVars = '%s\t%s';
formatSpecStep = 'Step Information: Rdamp=%f  (Run: %d/%d)';
formatSpecData = '%f\t%f';

    
varName = textscan(fileID, formatSpecVars, 1);  % pulls out variable name
stepName = textscan(fileID, formatSpecStep);    % pulls out run information 
stepR(1) = stepName{1,1};                           % pulls out run resistance
stepMax = stepName{1,3};                         % pulls out max step number
stepData = textscan(fileID, formatSpecData);
stepData = cell2mat(stepData);

% Time axis for interpolation of data
%timeRes = 50E-12; % 20GHz total BW (10Ghz positive)
%interpMethod = 'linear'; %pchip linear

% number of points (must round to real number
numIntPts = round(((stepData(end,1)-stepData(1,1))/timeRes)+1); % number of interpolated points

% time vector
interpTime = linspace(stepData(1,1), stepData(end,1), numIntPts);

% interpolates voltages
interpVolt = interp1(stepData(:,1),stepData(:,2),interpTime, interpMethod); %interpolates time

interpVoltAll = zeros(stepMax,numIntPts);
interpVoltAll(1, :) = interpVolt;

% loops through steps
if (1) %multiple steps
    for stepNum = 2:stepMax                             %start from 2 since run 1 above
        stepName = textscan(fileID, formatSpecStep);    % read step line
        stepR(stepNum) = stepName{1,1};                 % extracts resistance     
        stepData = textscan(fileID, formatSpecData);
        stepData = cell2mat(stepData);                  % converts to array

        interpVoltAll(stepNum, :) = interp1(stepData(:,1),stepData(:,2),interpTime);
   
    end
end

fclose(fileID);

end

