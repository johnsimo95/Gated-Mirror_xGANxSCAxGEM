  %% import comsol mirror TF 
  % The data that these call is in the GanLoadModelsData folder that I shared with you
 % the folder with the data and the folder with those scripts are linked to at the bottom of the simulation overview doc
filename='kimballCenterFieldTF_10MHz_10GHz.txt';
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


%%

% 
% %%
% highestFreq=5e10;

comsolFreqStep=fieldsVsFreq(2,1)-fieldsVsFreq(1,1);
posFreq=fieldsVsFreq(:,1);
% posFreqPad=transpose((posFreq(end)+comsolFreqStep):comsolFreqStep:highestFreq);
% posFreq=[posFreq;posFreqPad];
negFreq=-1*flip(posFreq);
mirrorTFpos=(fieldsVsFreq(:,5)./fieldsVsFreq(:,9));%/length(posFreq);
% mirrorTFposPad=mirrorTFpos(end)*ones(length(posFreqPad),1);
mirrorTFpos=[mirrorTFpos;mirrorTFposPad];
mirrorTFneg=flip(conj(mirrorTFpos));
comsolFreq=[posFreq;negFreq];
mirrorTF=[mirrorTFpos;mirrorTFneg];
% ALL ABOVE IS CONDITIONING RUN DATA FROM MATLAB pull into LIVESCRIPT
%%
%%%% Read in V over lumped element approx from SPICE - interpolate to same
%%%% frequencies as above^^^
% filenameSPICE='spiceExportRFmodelNoBearings1E8to1E10.txt';
% dataSPICE=readtable(filenameSPICE);
% dataSPICE=table2array(dataSPICE);
% time=dataSPICE(:,1);
% Vload=dataSPICE(:,2);
% timeReg=transpose(linspace(min(time),max(time),length(time))); %evenly spaced
% VloadReg=interp1(time,Vload,timeReg); %evenly spaced

workspaceName='200728 xGAN Symmetric Kimball Sweep.mat';
load(workspaceName);
r=7;%which damping resistor 
VloadReg=transpose(interpVoltAllDbl(r,:));
timeReg=transpose(interpTimeDbl);
% pading stuff to keep
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
% VloadRegPad=VloadReg;
% timeRegPad=timeReg;
Lpad=length(timeRegPad);

f = [0:Lpad/2-1, -Lpad/2:-1]*((1/(mean(diff(timeRegPad))))/Lpad);
%f = [0:Lpad/2-1, -Lpad/2:-1]*Fs; %??

fPos=f(1:length(f)/2);
dFreqSpice=mean(diff(fPos));

n=1;   %this up/downdamples things, don't need set n=1 keeps all the same
BW=1/T;
VloadTF=fft(VloadRegPad,length(f)*n);
% VloadTF=transpose(VloadTF);
freqStep=BW/(length(f)*n);

fPosDense=0:freqStep:BW/2-freqStep;
fNegDense=-BW/2:freqStep:0-freqStep;
fDense=transpose([fPosDense,fNegDense]);

VloadTFmatch=interp1(fDense,VloadTF,comsolFreq); % STOP USING THIS CODE HERE!! go to readLineDeflectionV4
mirrorTFmatch=interp1(comsolFreq,mirrorTF,fDense);

%check that FT of spice data works when you bring it back to time domain
Tcheck=1/range(comsolFreq);
IFFT_spiceFFT=ifft(VloadTFmatch,'symmetric');
%% freq response of whole system 
%%%%%%
% freqRespWholeSys=VloadTFmatch.*mirrorTF;
% T_wholeSys=1/range(comsolFreq);
freqRespWholeSys=VloadTF.*mirrorTFmatch;
T_wholeSys=1/range(fDense);

IFFTwholeSys=ifft(freqRespWholeSys,'symmetric');
fieldVsTime=IFFTwholeSys(padLength/2:(length(IFFTwholeSys)-padLength/2));

%%
set(0,'defaulttextInterpreter','latex') %latex axis labels

Legend=cell(length(stepR),1);
figure()
hold on
for r=1:length(stepR)
    Legend{r}=['R_{damp}=',num2str(stepR(r))];
    plot(interpTimeDbl,interpVoltAllDbl(r,:));
end
xlabel('time (s)');
ylabel('$V_{load}$');
legend(Legend);
grid on

figure()
plot(timeReg,VloadReg);
xlabel('time (s)');
ylabel('$V_{port}(t)$');
grid on

figure()
loglog(fDense,abs(VloadTF));
xlabel('Frequency (Hz)');
ylabel('$V_{port}(f)$');
grid on

figure()
loglog(fDense,abs(freqRespWholeSys));
xlabel('Frequency (Hz)')
ylabel('$V_{port}(f)\times TF_{mirror}(f)$');
grid on

figure()
plot(timeReg,fieldVsTime);
xlabel('time (s)');
ylabel('$|\vec{E}(x_0,y_0,z_0,t)|$');
grid on
%%
% figure
% plot(timeRegPad,VloadRegPad)
% hold on
% plot(0:Tcheck:Tcheck*(length(IFFT_spiceFFT)-1),IFFT_spiceFFT)
% xlabel('time');
% legend('TimeDomain Spice Signal','IFFT of FFT w/out >10GHz of time domain Spice')
% 
% 
% figure()
% semilogx(fDense,abs(VloadTF));
% xlabel('freq (Hz)')
% ylabel('V')
% title('Spice Switch Freq. Domain')
% 
% figure()
% % semilogx(comsolFreq,abs(freqRespWholeSys));
% semilogx(fDense,abs(freqRespWholeSys));
% 
% ylabel('Field Mag Response Whole Sys (V/m)')
% title('Transfer Function Whole System')
% hold on
% yyaxis right
% % semilogx(comsolFreq,unwrap(angle(freqRespWholeSys)));
% semilogx(fDense,unwrap(angle(freqRespWholeSys)));
% ylabel('Field Phase Response Whole Sys (rad)')
% xlabel('Hz')
% 
% figure()
% yyaxis left
% plot(0:T_wholeSys:T_wholeSys*(length(IFFTwholeSys)-1),IFFTwholeSys)
% xlabel('time (s)')
% ylabel('V/m at Mirror Center')
% hold on
% yyaxis right
% plot(timeRegPad,VloadRegPad);
% ylabel('V at Port)');

