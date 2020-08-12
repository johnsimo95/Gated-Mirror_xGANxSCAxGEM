% filename='portIV_freqResp_50ohm_realGeoNoBearings_wZ.txt';
filename='IV_GatedMirrorNoBearings.txt';

data=readtable(filename,'HeaderLines',5);
data=table2array(data);

freq=data(:,1);
freqSteps=diff(freq);
realV=data(:,2);
imagV=data(:,3)*1i;
realI=data(:,4);
imagI=data(:,5)*1i;
V=realV+imagV;
I=realI+imagI;
z=V./I;
%z=data(:,6);%V./I;
realZ=real(z);
imagZ=imag(z);
Zabs=abs(z);
ZPhase=unwrap(angle(z));

ratFit=rationalfit(freq,z);
sObject=sparameters(z,freq);
passiveRatFit=makepassive(ratFit,sObject);

% generateSPICE(ratFit,'testNetwork.ckt');

fitResp = freqresp(ratFit,freq);
fitRespMag=abs(fitResp);
fitRespPhase=angle(fitResp);

passiveTest=ispassive(ratFit);

passiveResp=freqresp(passiveRatFit,freq);
passiveRespMag=abs(passiveResp);
passiveRespPhase=angle(passiveResp);

figure()
semilogx(freq,Zabs);
hold on
semilogx(freq,fitRespMag);
% semilogx(freq,passiveRespMag);
hold on
title('Magnitude');

figure()
semilogx(freq,ZPhase);
hold on
semilogx(freq,fitRespPhase);
% semilogx(freq,passiveRespPhase);
title('Phase');
% %%
% data=frd(z,freq);
% TFsysTest=tfest(data,7);
% 
% zTF=tf(TFsysTest.Numerator,TFsysTest.Denominator);
% %bode(zTF);
% fitTFatFreq=freqresp(zTF,freq);
% fitTFatFreq=squeeze(fitTFatFreq);
% magFitZ=abs(fitTFatFreq);
% phaseFitZ=imag(fitTFatFreq);
% 
% figure()
% title('Magnitude')
% semilogx(freq,Zabs)
% hold on
% yyaxis right
% semilogx(freq,magFitZ);
% 
% figure()
% title('phase')
% semilogx(freq,ZPhase)
% hold on
% yyaxis right
% semilogx(freq,phaseFitZ);

