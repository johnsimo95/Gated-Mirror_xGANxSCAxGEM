addpath('RawData');
filename='portIV_freqResp_50ohm_realGeoNoBearins.txt';

data=readtable(filename,'HeaderLines',5);
data=table2array(data);

freq=data(:,1);
realV=data(:,2);
imagV=data(:,3)*1i;
realI=data(:,4);
imagI=data(:,5)*1i;
V=realV+imagV;
I=realI+imagI;
z=V./I;
realZ=real(z);
imagZ=imag(z);
Zabs=abs(z);
Zphase=(angle(z));

figure()
plot(freq,Zabs)
figure()
plot(freq,Zphase)

%%
order=11;%floor((length(freq)-1)/2);

%see if it works with just linearly scaling frequency axis
freqScale=freq.*1E-9;
scaleRealZP=polyfit(freqScale,realZ,order);
scaleRealZFit=polyval(scaleRealZP,freqScale);
scaleImagZP=polyfit(freqScale,imagZ,order);
scaleImagZFit=polyval(scaleImagZP,freqScale);

figure(1)
hold off
plot(freqScale,imagZ)
hold on
yyaxis right
plot(freqScale,scaleImagZFit)
hold off
%%
%translating frequency axis into units of stddev


[realZP_trans,realS,realMu]=polyfit(freq,realZ,order);
[imagZP_trans,imagS,imagMu]=polyfit(freq,imagZ,order);

realZFit_trans=polyval(realZP_trans,freq,realS,realMu);
imagZFit_trans=polyval(imagZP_trans,freq,imagS,imagMu);

figure(3)
hold off
plot(freq,realZ)
hold on
yyaxis right
plot(freq,realZFit_trans);
hold off

% normRealZ=1/realMu(2)-realMu(1)/realMu(2);
% what_RealZ=sym2poly(subs(poly2sym(realP),poly2sym(normRealZ)));
% coeffs_freqDomain=polyval(what_RealZ,freq);
