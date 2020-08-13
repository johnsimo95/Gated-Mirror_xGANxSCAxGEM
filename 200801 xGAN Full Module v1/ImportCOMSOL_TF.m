function [fieldsVsFreq, mirrorTF, comsolFreq] = ImportCOMSOL_TF(filename, makeSymmetric)
%IMPORTCOMSOL_TF Summary of this function goes here
%   Detailed explanation goes here

    fileID=fopen(filename);
    headers=textscan(fileID,'%s',1,'delimiter','\n','HeaderLines',8);
    headers=char(headers{1});
    fseek(fileID,0,'bof');
    rawPoints=textscan(fileID,'%f','HeaderLines',9);
    rawPoints=rawPoints{1};
    rawPoints=rawPoints(4:end);
    fclose(fileID);
    freqIndex=strfind(headers,'=');
    freqSteps=zeros(length(freqIndex),1);%sscanf(headers(timeIndex(1)+1:end),'%f');
    for n=1:length(freqIndex)
        if freqIndex(n)<(length(headers)-20)
            freqSteps(n)=sscanf(headers(freqIndex(n)+1:freqIndex(n)+20),'%f',1);
        else
            freqSteps(n)=sscanf(headers(freqIndex(n)+1:end),'%f',1);
        end
    end
    fieldsVsFreq=zeros(length(rawPoints)/7,9); %time,Ex,Ey,Ez,|E|,Vport,Vnode
    perpPlaneFieldRF=zeros(length(rawPoints)/7,1); %magnitude of field components in plane perpendicular to optical axis (y-z)
    n=1;
    k=1;
    while k<=(length(rawPoints)/7)
        fieldsVsFreq(k,1)=freqSteps(n);
        fieldsVsFreq(k,2)=rawPoints(n); %ex
        fieldsVsFreq(k,3)=rawPoints(n+1); %ey
        fieldsVsFreq(k,4)=rawPoints(n+2); %ez
        fieldsVsFreq(k,5)=sqrt(fieldsVsFreq(k,2)^2+fieldsVsFreq(k,3)^2+fieldsVsFreq(k,4)^2); %e magnitude
        fieldsVsFreq(k,6)=rawPoints(n+3); %node 1 voltage
        fieldsVsFreq(k,7)=rawPoints(n+4); %ideal voltage supply
        fieldsVsFreq(k,8)=rawPoints(n+5); %node 2 voltage
        fieldsVsFreq(k,9)=rawPoints(n+6); %voltage at port
        perpPlaneFieldRF(k)=sqrt(abs(fieldsVsFreq(k,3))^2+abs(fieldsVsFreq(k,4))^2);
        n=n+7; 
        k=k+1;
    end

    %Default state of mirror!
    mirrorTF = 0;
    comsolFreq = 0;
    if (makeSymmetric)
        comsolFreqStep=fieldsVsFreq(2,1)-fieldsVsFreq(1,1)
        posFreq=fieldsVsFreq(:,1);
        % posFreqPad=transpose((posFreq(end)+comsolFreqStep):comsolFreqStep:highestFreq);
        % posFreq=[posFreq;posFreqPad];
        negFreq=-1*flip(posFreq);
        mirrorTFpos=(fieldsVsFreq(:,5)./fieldsVsFreq(:,9));%/length(posFreq);
        % mirrorTFposPad=mirrorTFpos(end)*ones(length(posFreqPad),1);
        mirrorTFpos=[mirrorTFpos];
        mirrorTFneg=flip(conj(mirrorTFpos));
        comsolFreq=[posFreq;negFreq];
        mirrorTF=[mirrorTFpos;mirrorTFneg];
    end

end

