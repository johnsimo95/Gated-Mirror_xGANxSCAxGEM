%constants/params
k=1e3;
q=1.6e-19; %charge of electron.
eVtoJ=q;
m= 9.1e-31; %mass of electron
E0=10*k*eVtoJ;
vel=sqrt(2*E0/m); %velocity for force due to magnetic field calcs

filename='deflectionLineTest.txt';
fileID=fopen(filename);
headers=textscan(fileID,'%s',1,'delimiter','\n','HeaderLines',8);
headers=char(headers{1});
fseek(fileID,0,'bof'); %point back to start of file
rawPoints=textscan(fileID,'%f','HeaderLines',9);
rawPoints=rawPoints{1};
% rawPoints=rawPoints(4:end);
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
freqSteps=unique(freqSteps);
negFreq=-1*flip(freqSteps);
posFreq=freqSteps-mean(diff(freqSteps));
comsolFreq=[freqSteps;negFreq];
comsolFreq=unique(comsolFreq);
numPosSteps=length(rawPoints)/(5*length(freqSteps)+3);
freqDomainFieldTFsOfX=zeros(length(freqSteps),4,numPosSteps);
%Ey/V for all x, Ez/V for allx, By/V for all x, Bz/V for allx at each
%freqency
%
forceOverVin_vsFreq=zeros(length(freqSteps),4); %integral of Ey,Ez,By,Bz dz for each freq (*q or *qv)
%^ result should be F(f)/Vin(f)

for n=1:length(freqSteps)

      temp=zeros(numPosSteps,6); %X,Ey,Ez,By,Bz,Vport
  for k=1:numPosSteps
      temp(k,1)=rawPoints(1+(k-1)*(5*length(freqSteps)+3)); %X
      temp(k,2)=rawPoints(1+3+(n-1)*5+(k-1)*(5*length(freqSteps)+3)); %Ey
      temp(k,3)=rawPoints(1+4+(n-1)*5+(k-1)*(5*length(freqSteps)+3)); %Ez
      temp(k,4)=rawPoints(1+5+(n-1)*5+(k-1)*(5*length(freqSteps)+3)); %By
      temp(k,5)=rawPoints(1+6+(n-1)*5+(k-1)*(5*length(freqSteps)+3)); %Bz
      temp(k,6)=rawPoints(1+7+(n-1)*5+(k-1)*(5*length(freqSteps)+3)); %Vport
  end

  pos=temp(:,1);
  %need to integrate over VxB*dx and E/q*dx for each freq 
  % actually get TF of these compared to port voltage?
  EyOverV=temp(:,2)./temp(:,6);
  EzOverV=temp(:,3)./temp(:,6);
  ByOverV=temp(:,4)./temp(:,6);
  BzOverV=temp(:,5)./temp(:,6);
 
  freqDomainFieldTFsOfX(n,1,:)=EyOverV;
  freqDomainFieldTFsOfX(n,2,:)=EzOverV;
  freqDomainFieldTFsOfX(n,3,:)=ByOverV;
  freqDomainFieldTFsOfX(n,4,:)=BzOverV;

%   forceOverVin_vsFreq(n,1)=trapz(pos,EyOverV*q);
%   forceOverVin_vsFreq(n,2)=trapz(pos,EzOverV*q);
%   forceOverVin_vsFreq(n,3)=trapz(pos,ByOverV*q*vel); %to get force in y, shouldn't I do Bz, since cross product?
%   forceOverVin_vsFreq(n,4)=trapz(pos,BzOverV*q*vel); %^^^
  
%   elecForce=trapz(temp(:,4)./temp(:,8)*q);
%   magForce=trapz(temp(:,7)./temp(:,8)*vel*q);
%   sortedData{n,3}=[elecForce,magForce];
%   sortedData{n,4}=elecForce+magForce;
end

freqDomainFieldTFsOfX_negFreq=zeros(length(freqSteps),4,numPosSteps); %get the tf of each field component at each point in space to have the negative frequency part
fieldTFperX=zeros(2*length(freqSteps),4,numPosSteps); %for combining into full double sided TF

for p=1:numPosSteps
   freqDomainFieldTFsOfX_negFreq(:,1,p)=flip(conj(freqDomainFieldTFsOfX(:,1,p)));
   freqDomainFieldTFsOfX_negFreq(:,2,p)=flip(conj(freqDomainFieldTFsOfX(:,2,p)));
   freqDomainFieldTFsOfX_negFreq(:,3,p)=flip(conj(freqDomainFieldTFsOfX(:,3,p)));
   freqDomainFieldTFsOfX_negFreq(:,4,p)=flip(conj(freqDomainFieldTFsOfX(:,4,p)));
   
   fieldTFperX(:,1,p)=[freqDomainFieldTFsOfX(:,1,p);freqDomainFieldTFsOfX_negFreq(:,1,p)];
   fieldTFperX(:,2,p)=[freqDomainFieldTFsOfX(:,2,p);freqDomainFieldTFsOfX_negFreq(:,2,p)];
   fieldTFperX(:,3,p)=[freqDomainFieldTFsOfX(:,3,p);freqDomainFieldTFsOfX_negFreq(:,3,p)];
   fieldTFperX(:,4,p)=[freqDomainFieldTFsOfX(:,4,p);freqDomainFieldTFsOfX_negFreq(:,4,p)];
   
end


%% Can take all of this! Import all below from readline defclection Take 94-119
workspaceName='200728 xGAN Symmetric Kimball Sweep.mat';
load(workspaceName);
r=7;%which damping resistor 
vLoadRaw=transpose(interpVoltAllDbl(r,:));
timeRaw=transpose(interpTimeDbl);
% time=downsample(timeRaw,5);
% vLoad=downsample(vLoadRaw,5);
time=decimate(timeRaw,7);
vLoad=decimate(vLoadRaw,7);

T=mean(diff(time)); %time iterations
BW=1/(2*T);
Fs=1/T; %sampling frequency
L=length(time);
padLength=floor(L/2)*100;
[trash,midVidx]=min(abs(vLoad-mean(vLoad)));
shift=floor(L/2)-midVidx;
padFront=vLoad(1)*ones(shift+padLength/2,1);
padEnd=vLoad(end)*ones(-shift+padLength/2,1);
vLoadPad=[padFront;vLoad;padEnd];
timePad=T*transpose(0:length(vLoadPad)-1);
Lpad=length(timePad);
fPos=(0:(Lpad/2-1))*(Fs/Lpad);
fNeg=((-Lpad/2):-1)*(Fs/Lpad);
f=transpose([0:Lpad/2-1, -Lpad/2:-1]*(Fs/Lpad));
vLoadTF=fft(vLoadPad,length(f));
%%
%find actual response to field
if max(f)> max(comsolFreq)
    vLoadTFmatch=interp1(f,vLoadTF,comsolFreq); 
    vLoadMatchCheck=ifft(vLoadTFmatch,'symmetric');

    fieldRespToSpice=zeros(2*length(freqSteps),4,numPosSteps); %for combining into full double sided TF

    for p=1:numPosSteps
        fieldRespToSpice(:,1,p)=fieldTFperX(:,1,p).*vLoadTFmatch;
        fieldRespToSpice(:,2,p)=fieldTFperX(:,2,p).*vLoadTFmatch;
        fieldRespToSpice(:,3,p)=fieldTFperX(:,3,p).*vLoadTFmatch;
        fieldRespToSpice(:,4,p)=fieldTFperX(:,4,p).*vLoadTFmatch;
    end

    fieldRespToSpice_posFreq=fieldRespToSpice(1:length(freqSteps),:,:);
    % forceOfFreq=vLoadTFmatch.*forceOverVin_vsFreq;
    % forceOfTime=ifft(forceOfFreq,'symmetric');
    % T_forForce=1/range(freqSteps);
    % timeForForce=[0:T_forForce:T_forForce*(length(forceOfTime)-1)];
    % vLoadTF_IFFT=ifft(vLoadTFmatch,'symmetric');

    T_inc=1/range(comsolFreq);
    t=[0:T_inc:T_inc*(length(comsolFreq)-1)];
    deltaX=mean(diff(pos));
    deltaW=2*pi*mean(diff(freqSteps));

else
    fieldRespToSpice_matchCtoS=zeros(length(f),4,numPosSteps); %interpolation stuff, tricky!
    fieldTFperX_match=zeros(length(f),4,numPosSteps);
    for p=1:numPosSteps
        fieldTFperX_match(:,1,p)=interp1(comsolFreq,fieldTFperX(:,1,p),f);
        fieldTFperX_match(:,2,p)=interp1(comsolFreq,fieldTFperX(:,2,p),f);
        fieldTFperX_match(:,3,p)=interp1(comsolFreq,fieldTFperX(:,3,p),f);
        fieldTFperX_match(:,4,p)=interp1(comsolFreq,fieldTFperX(:,4,p),f);

        
        fieldRespToSpice_matchCtoS(:,1,p)=fieldTFperX_match(:,1,p).*vLoadTF;
        fieldRespToSpice_matchCtoS(:,2,p)=fieldTFperX_match(:,2,p).*vLoadTF;
        fieldRespToSpice_matchCtoS(:,3,p)=fieldTFperX_match(:,3,p).*vLoadTF;
        fieldRespToSpice_matchCtoS(:,4,p)=fieldTFperX_match(:,4,p).*vLoadTF;
    end
    T_inc=1/range(f);
    t=[0:T_inc:T_inc*(length(f)-1)];
    deltaX=mean(diff(pos));
    deltaW=2*pi*mean(diff(fPos));

    fieldRespToSpice_posFreq=fieldRespToSpice_matchCtoS(1:length(fPos),:,:);

end

%%
% syms t
if max(f)>max(comsolFreq)
    %need to replace freqDomainFieldsTFofX with that*FFT of spicedomain Signal
    intFdx_perFreq=zeros(length(freqSteps),4,length(t)); %(sum of F(x,w)*dx over all x)*exp(jwt), at each w, how to get t???
    %intFdx_perFreq=sym('F',[length(freqSteps) 4]);
    for n=1:length(freqSteps)
        intFdx_perFreq(n,1,:)=sum(fieldRespToSpice_posFreq(n,1,:))*deltaX*exp(j*freqSteps(n)*2*pi*t);
        intFdx_perFreq(n,2,:)=sum(fieldRespToSpice_posFreq(n,2,:))*deltaX*exp(j*freqSteps(n)*2*pi*t);
        intFdx_perFreq(n,3,:)=sum(fieldRespToSpice_posFreq(n,3,:))*deltaX*exp(j*freqSteps(n)*2*pi*t);  
        intFdx_perFreq(n,4,:)=sum(fieldRespToSpice_posFreq(n,4,:))*deltaX*exp(j*freqSteps(n)*2*pi*t);
    end

    intFdx=zeros(length(t),4);
    %intFdx=sym('Fdx',[1 4]);
    for tstep=1:length(t)
       intFdx(tstep,1)=sum(intFdx_perFreq(:,1,tstep));
       intFdx(tstep,2)=sum(intFdx_perFreq(:,2,tstep));
       intFdx(tstep,3)=sum(intFdx_perFreq(:,3,tstep));   
       intFdx(tstep,4)=sum(intFdx_perFreq(:,4,tstep));
    end

    % intFdx(:,1)=sum(intFdx_perFreq(:,1,:))*deltaW;
    % intFdx(:,2)=sum(intFdx_perFreq(:,2,:))*deltaW;
    % intFdx(:,3)=sum(intFdx_perFreq(:,3,:))*deltaW;
    % intFdx(:,4)=sum(intFdx_perFreq(:,4,:))*deltaW;
    % 
else
        intFdx_perFreq=zeros(length(fPos),4,length(t)); %(sum of F(x,w)*dx over all x)*exp(jwt), at each w, how to get t???
    %intFdx_perFreq=sym('F',[length(freqSteps) 4]);
    for n=1:length(fPos)
        intFdx_perFreq(n,1,:)=sum(fieldRespToSpice_posFreq(n,1,:))*deltaX*exp(j*fPos(n)*2*pi*t);
        intFdx_perFreq(n,2,:)=sum(fieldRespToSpice_posFreq(n,2,:))*deltaX*exp(j*fPos(n)*2*pi*t);
        intFdx_perFreq(n,3,:)=sum(fieldRespToSpice_posFreq(n,3,:))*deltaX*exp(j*fPos(n)*2*pi*t);  
        intFdx_perFreq(n,4,:)=sum(fieldRespToSpice_posFreq(n,4,:))*deltaX*exp(j*fPos(n)*2*pi*t);
    end

    intFdx=zeros(length(t),4);
    %intFdx=sym('Fdx',[1 4]);
    for tstep=1:length(t)
       intFdx(tstep,1)=sum(intFdx_perFreq(:,1,tstep));
       intFdx(tstep,2)=sum(intFdx_perFreq(:,2,tstep));
       intFdx(tstep,3)=sum(intFdx_perFreq(:,3,tstep));   
       intFdx(tstep,4)=sum(intFdx_perFreq(:,4,tstep));
    end
    
end
%%
set(0,'defaulttextInterpreter','latex') %latex axis labels
if max(f)>max(comsolFreq)
figure
plot(time,vLoad)
hold on
plot(t,abs(vLoadMatchCheck));
title('$V_{port}(t)$')
legend('V(t)','iFFT of FFT of downsampled V(t)');
figure
plot(f,abs(vLoadTF))
hold on
plot(comsolFreq,abs(vLoadTFmatch))
legend('$V_{port}(f)$','$FFT of downsampled V_{port}(t)');

else
    figure
    plot(timePad,vLoadPad);
    title('$V_{port}(t)$');
end

figure
plot(t,abs(intFdx(:,1)));
title('$F_{Ey}(t)$');
figure
plot(t,abs(intFdx(:,2)));
title('$F_{Ez}(t)$');
figure
plot(t,abs(intFdx(:,3)));
title('$F_{By}(t)$');
figure
plot(t,abs(intFdx(:,4)));
title('$F_{Bz}(t)$');

%%
figure
plot(time,abs(intFdx(padLength/2:end-padLength/2-1,1)));
title('$F_{Ey}(t)$');
figure
plot(time,abs(intFdx(padLength/2:end-padLength/2-1,2)));
title('$F_{Ez}(t)$');
figure
plot(time,abs(intFdx(padLength/2:end-padLength/2-1,3)));
title('$F_{By}(t)$');
figure
plot(time,abs(intFdx(padLength/2:end-padLength/2-1,4)));
title('$F_{Bz}(t)$');

