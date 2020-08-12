%% high frequency
filename='kimballMirrorFieldsTF_R50_HiFreq.txt';
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
%% low frequency
filename='kimballMirrorFieldsTF_R50.txt';
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
fieldsVsFreq=[fieldsVsFreqLF;fieldsVsFreqRF];
perpPlaneField=[perpPlaneFieldLF;perpPlaneFieldRF];

fieldsVsFreqMag=abs(fieldsVsFreq(:,5));
sourceVsFreqMag=abs(fieldsVsFreq(:,7));
mirrorTF=fieldsVsFreq(:,5)./fieldsVsFreq(:,7);
mirrorTFmag=abs(mirrorTF);
mirrorTFmagDB=20*log10(mirrorTFmag);

figure()
yyaxis left
semilogx(fieldsVsFreq(:,1),fieldsVsFreqMag,'x');
xlabel('freq (Hz)');
ylabel('Field Strength (V/m)');
yyaxis right
hold on
semilogx(fieldsVsFreq(:,1),mirrorTFmagDB);
ylabel('Mirror TF (dB per m)');
grid on

figure()
semilogx(fieldsVsFreq(:,1),perpPlaneField,'x')
xlabel('freq (Hz)');
ylabel('Field Strength in Perp. Plane (V/m)');
grid on
