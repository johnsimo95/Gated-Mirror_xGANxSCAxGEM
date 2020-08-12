


filenameSPICE='SampleDataExport.txt';
dataSPICE=readtable(filenameSPICE);
dataSPICE=table2array(dataSPICE);
time=dataSPICE(:,1);
Vload=dataSPICE(:,2);
timeReg=linspace(min(time),max(time),length(time)); %evenly spaced
VloadReg=interp1(time,Vload,timeReg); %evenly spaced
T=mean(diff(timeReg)); %time iterations
Fs=1/T; %sampling frequency
L=length(timeReg);
f=Fs*(0:(L/2))/L;
VloadTF=fft(Vload);
P2=abs(VloadTF/L);
P1=P2(1:L/2+1);
P1(2:end-1)=2*P1(2:end-1);

figure()
semilogx(f,P1)
xlabel('f (Hz)');
ylabel('Vload');