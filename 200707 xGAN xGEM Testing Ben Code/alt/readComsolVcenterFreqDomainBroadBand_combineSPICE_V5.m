%% RF 
filename='RFfreqResp_50ohm_NoBearingsCenter.txt';
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
fieldsVsFreqRF=fieldsVsFreq;
%% LF
filename='LFfreqResp_NoBearings50ohmCenter_vFine_widerBoundary.txt';
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
perpPlaneFieldLF=zeros(length(rawPoints)/7,1); %magnitude of field components in plane perpendicular to optical axis (y-z)
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
    perpPlaneFieldLF(k)=sqrt(abs(fieldsVsFreq(k,3))^2+abs(fieldsVsFreq(k,4))^2);
    n=n+7; 
    k=k+1;
end
fieldsVsFreqLF=fieldsVsFreq;
%%
% RF just to get the high end of the spice response to match up with Comsol
% TF
filename='RFfreqResp_NoBearings50OhmCenter_HFspice.txt';
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
perpPlaneFieldRF_HFspice=zeros(length(rawPoints)/7,1); %magnitude of field components in plane perpendicular to optical axis (y-z)
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
    perpPlaneFieldRF_HFspice(k)=sqrt(abs(fieldsVsFreq(k,3))^2+abs(fieldsVsFreq(k,4))^2);
    n=n+7; 
    k=k+1;
end
fieldsVsFreqRF_HFspice=fieldsVsFreq;

%%
% scrap together all the different Field vs Port V transfer functions:
% Low Freq + RF + RF adjusted for the high end of the spice sim
if fieldsVsFreqLF(end,1)==fieldsVsFreqRF(1,1)
    fieldsVsFreq=[fieldsVsFreqLF(1:end-1,:);fieldsVsFreqRF;fieldsVsFreqRF_HFspice(2:end,:)];
    perpPlaneField=[perpPlaneFieldLF(1:end-1,:);perpPlaneFieldRF;perpPlaneFieldRF_HFspice(2:end,:)];
else
    fieldsVsFreq=[fieldsVsFreqLF;fieldsVsFreqRF;fieldsVsFreqRF_HFspice(2:end,:)];
    perpPlaneField=[perpPlaneFieldLF;perpPlaneFieldRF;perpPlaneFieldRF_HFspice(2:end,:)];
end
posFreq=fieldsVsFreq(:,1);
negFreq=-1*flip(posFreq);
mirrorTFpos=(fieldsVsFreq(:,5)./fieldsVsFreq(:,7));%/length(posFreq);
mirrorTFneg=flip(conj(mirrorTFpos));
comsolFreq=[posFreq;negFreq];
mirrorTF=[mirrorTFpos;mirrorTFneg];

%%
%%%% Read in V over lumped element approx from SPICE - interpolate to same
%%%% frequencies as above^^^
filenameSPICE='spiceExportRFmodelNoBearings1E8to1E10.txt';
dataSPICE=readtable(filenameSPICE);
dataSPICE=table2array(dataSPICE);
time=dataSPICE(:,1);
Vload=dataSPICE(:,2);
timeReg=transpose(linspace(min(time),max(time),length(time))); %evenly spaced
VloadReg=interp1(time,Vload,timeReg); %evenly spaced
T=mean(diff(timeReg)); %time iterations
Fs=1/T; %sampling frequency
L=length(timeReg);

padLength=floor(L/2)*100;
[trash,midVidx]=min(abs(VloadReg-mean(VloadReg)));
shift=floor(L/2)-midVidx;
padFront=VloadReg(1)*ones(shift+padLength/2,1);
padEnd=VloadReg(end)*ones(-shift+padLength/2,1);

VloadRegPad=[padFront;VloadReg;padEnd];
timeRegPad=T*transpose(0:length(VloadRegPad)-1);
Lpad=length(timeRegPad);

f = [0:Lpad/2-1, -Lpad/2:-1]*((1/(mean(diff(timeRegPad))))/Lpad);
fPos=f(1:length(f)/2);
dFreqSpice=mean(diff(fPos));

n=1;
BW=1/mean(diff(timeRegPad));
VloadTF=fft(VloadRegPad,length(f)*n);
% VloadTF=transpose(VloadTF);
freqStep=BW/(length(f)*n);

fPosDense=0:freqStep:BW/2-freqStep;
fNegDense=-BW/2:freqStep:0-freqStep;
fDense=transpose([fPosDense,fNegDense]);

mirrorTFmatch=interp1(comsolFreq,mirrorTF,fDense);
%% freq response of whole system and AWG model
%%%%%%
lDense=length(fDense);
freqRespWholeSys=VloadTF.*mirrorTFmatch;

invMirrorTFmatch=mirrorTFmatch.^-1;
iMTFmData=frd(invMirrorTFmatch(1:lDense/2),fDense(1:lDense/2),'FrequencyUnit','Hz');
invMirrorFit=tfest(iMTFmData,15);
invMirrorFitMatch=squeeze(freqresp(invMirrorFit,fDense(1:lDense/2),'Hz'));

%%
figure()
semilogx(comsolFreq,mirrorTF);
hold on
semilogx(fDense,mirrorTFmatch);
xlabel('Freq (Hz)')
ylabel('(V/m)/(V_{port})')
legend('Mirror TF from Comsol','Mirror TF interpolated')

figure()
semilogx(fDense,abs(VloadTF));
xlabel('freq (Hz)')
ylabel('V')
title('Spice Switch Freq. Domain')

figure()
semilogx(fDense,abs(freqRespWholeSys));
ylabel('Field Mag Response Whole Sys (V/m)')
title('Transfer Function Whole System')
hold on
yyaxis right
semilogx(fDense,unwrap(angle(freqRespWholeSys)));
ylabel('Field Phase Response Whole Sys (rad)')
xlabel('Hz')


T=1/freqStep;%(mean(diff(posFreqLin)));
T=(1/(range(fDense)));%/length(fDense);
brokenTF=find(isnan(freqRespWholeSys)==1);
freqRespWholeSys(brokenTF)=[];
IFFTwholeSys=ifft(freqRespWholeSys,'symmetric');
figure()
plot(0:T:T*(length(IFFTwholeSys)-1),IFFTwholeSys)
xlabel('time (s)')
ylabel('V/m')
title('E-Field at Mirror Center');
