function [trFnRaw, trFnSymmetric] = ImportCOMSOL_TF(filename, portVoltage, makeSymmetric)
%IMPORTCOMSOL_TF Summary of this function goes here
%   Detailed explanation goes here

    data=readtable(filename,'HeaderLines',5);
    data=table2array(data);

    trFnRaw.freq = data(:,1);
    trFnRaw.fieldTF = data(:,2)./portVoltage;
    %Default state of mirror!
    
    if (makeSymmetric)
        comsolFreqStep=trFnRaw.freq(2,1)-trFnRaw.freq(1,1)
        posFreq=trFnRaw.freq;
        negFreq=-1*flip(posFreq);
        trFnSymmetric.freq = [posFreq;negFreq];

        mirrorTFpos=(trFnRaw.fieldTF);%/length(posFreq);
        mirrorTFneg=flip(conj(mirrorTFpos));
        trFnSymmetric.fieldTF = [mirrorTFpos;mirrorTFneg];
    end

end

