% filename='portIV_freqResp_50ohm_realGeoNoBearings_wZ.txt';
% filename='IV_GatedMirrorNoBearings.txt';



%% EXTRACTING AND COMBINING RAW DATA
%% from high Frequency
filenameRF='portIVkimballMirror_R50_HiFreq.txt';
dataRF=readtable(filenameRF,'HeaderLines',5);
dataRF=table2array(dataRF);
freqRF=dataRF(:,1);
% freqStepsRF=diffRF(freq);
realVRF=dataRF(:,2);
imagVRF=dataRF(:,3)*1i;
realIRF=dataRF(:,4);
imagIRF=dataRF(:,5)*1i;
VRF=realVRF+imagVRF;
IRF=realIRF+imagIRF;
zRF=VRF./IRF;


%% From low Frequency
filenameLF='portIVkimballMirror_R50.txt';
dataLF=readtable(filenameLF,'HeaderLines',5);
dataLF=table2array(dataLF);
freqLF=dataLF(:,1);
% freqStepsLF=diffLF(freq);
realVLF=dataLF(:,2);
imagVLF=dataLF(:,3)*1i;
realILF=dataLF(:,4);
imagILF=dataLF(:,5)*1i;
VLF=realVLF+imagVLF;
ILF=realILF+imagILF;
zLF=VLF./ILF;


%% Combine data
if dataLF(end,1)==dataRF(1,1)
    data=[dataLF(1:end-1,:);dataRF];
else
    data=[dataLF;dataRF];
end
freq=data(:,1);
V=data(:,2)+data(:,3)*1j;
I=data(:,4)+data(:,5)*1j;
z=V./I;



%%
%    z=zRF;
%    freq=freqRF;  
%%
Zabs=abs(z);
ZPhase=unwrap(angle(z));
z0=50; %resistance in source
% Vs=data(:,9);
% Pavs=(Vs.^2)/(8*z0);
% Pin=0.5*(V.*conj(I)+conj(V).*I);
% Pr=Pavs-Pin;
% Gamma=sqrt(Pr./Pavs);
Gamma=z2gamma(z);
    % sObj=sparameters(Gamma,freq);
    % ratFit=rationalfit(sObj);
Gamma=makepassive(Gamma);
ratFit=rationalfit(freq,squeeze(Gamma));%,'TendsToZero',true);
%%
%generateSPICE(ratFit,'inductor_10MHzto10GHz.ckt');
%%
fitResp = freqresp(ratFit,freq);
fitRespMag=abs(fitResp);
fitRespPhase=angle(fitResp);
%%
figure()
yyaxis left
semilogx(freq,Zabs);
ylabel('Magnitude (\Omega)');
xlabel('Hz')
hold on
yyaxis right
semilogx(freq,ZPhase);
ylabel('Phase (rad)');
title('Effective Impedance at Port');
%%
figure()
yyaxis left
semilogx(freq,Zabs);
ylabel('Simulated impedance magnitude');
hold on
yyaxis right
semilogx(freq,fitRespMag);
% semilogx(freq,passiveRespMag);
title('Magnitude');
xlabel('Hz')
ylabel('S-param rat. fit Impedance Magnitude');
% 
% figure()
% semilogx(freq,ZPhase);
% ylabel('Simulated impedance phase');
% hold on
% yyaxis right
% semilogx(freq,fitRespPhase);
% % semilogx(freq,passiveRespPhase);
% title('Phase');
% xlabel('Hz')
% ylabel('S-param rat. fit impedance phase');
