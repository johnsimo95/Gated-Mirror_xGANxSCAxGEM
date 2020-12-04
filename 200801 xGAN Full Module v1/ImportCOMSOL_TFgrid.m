function [trFnRawGrid, trFnSymmetricGrid] = ImportCOMSOL_TFgrid(filename, portVoltage, makeSymmetric)
%IMPORTCOMSOL_TF Summary of this function goes here
%   Detailed explanation goes here

    data=readtable(filename,'HeaderLines',9);
    data=table2array(data);
    sizeArr = size(data);
    sizeArr = sizeArr(1,1);
    % NO DC TERM IN RAW!
    for(ii = 1:numPts)
        freq = freq.min:freq.step:freq.max;
        trFnRaw(ii).x = (data(:,1)).';
        trFnRaw(ii).y = (data(:,2)).';
        trFnRaw(ii).z = (data(:,3)).';
        trFnRaw(ii).fieldTF = (data(:,4:end)./portVoltage).';
    end
    %Default state of mirror!
    
    if (makeSymmetric)
        comsolFreqStep=trFnRaw.freq(1,2)-trFnRaw.freq(1,1);
        posFreq=trFnRaw.freq;
        negFreq=-1*flip(posFreq);
        % Ordering of FFT is 0-Pos-Neg
        trFnSymmetric.freq = [0,posFreq,negFreq]; %add DC term!

        mirrorTFpos=(trFnRaw.fieldTF);%/length(posFreq);
        mirrorTFneg=flip(conj(mirrorTFpos)); %new side conjugate symmetric!!
        trFnSymmetric.fieldTF = [0, mirrorTFpos, mirrorTFneg]; %setting DC term to 0! AC coupling...
    end

end

