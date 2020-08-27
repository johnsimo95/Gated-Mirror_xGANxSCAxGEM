function [trFnRaw, trFnSymmetric] = ImportCOMSOL_TF(filename, portVoltage, makeSymmetric)
%IMPORTCOMSOL_TF Summary of this function goes here
%   Detailed explanation goes here

    data=readtable(filename,'HeaderLines',5);
    data=table2array(data);

    % NO DC TERM IN RAW!
    trFnRaw.freq = (data(:,1)).';
    trFnRaw.fieldTF = (data(:,2)./portVoltage).';
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

