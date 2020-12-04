function [ sinParams ] = streakImgRead( fileID )

imgArray=imread(fileID);
imgArray=rgb2gray(imgArray);

xIndex=transpose(1:length(imgArray(1,:)));
xLength=length(xIndex);
%at this point, we would multiply xLength by the calibration factor for
%distance/convert from distance to time
yIndex=1:length(imgArray(:,1));
yLength=length(yIndex);

%%

traceYindex=zeros(xLength,1);
maxYindex=cell(xLength,1);

for n=1:xLength
   ySlice=imgArray(:,n);
   maxY=max(ySlice);
   maxYindex{n}=find(ySlice==maxY);
   
%    if length(maxYindex)==1
%        traceYindex(n)=maxYindex(n);
%    else
       traceYindex(n)=mean(maxYindex{n});
%    end
    
end
%at this point we would multiply traceYindex by the calibration factor for
%Field strength per pixel (index and pixel are interchangeable here)
%%
%attempt to clean up weird spikes 
yMaxNormal=fitdist(traceYindex,'Normal');

for n=1:length(traceYindex)
   if (traceYindex(n)<yMaxNormal.mu-yMaxNormal.sigma)||(traceYindex(n)>yMaxNormal.mu+yMaxNormal.sigma)
    traceYindex(n)=mean([traceYindex(n-1),traceYindex(n+1)]);
   end
end
traceYindexSmoothed=smooth(traceYindex,0.02); %play around with smoothing percentage (second number)

%going to try to fit only the middle fifty percent of points b/c it looks
%shitty at the edges of the image
lineFit=polyfit(xIndex(floor(length(xIndex)/4):length(xIndex)-floor(length(xIndex)/4)),traceYindexSmoothed(floor(length(xIndex)/4):length(xIndex)-floor(length(xIndex)/4)),1);


slope=lineFit(1);
yInt=lineFit(2);
yFit=yInt+slope*xIndex;
traceYindexSmoothedCentered=traceYindexSmoothed-yInt;
theta=-atand(slope);
%%
%coordinate transformation so data doesn't look like x*sin(x), but just
%sin(x) centered around y=0
yIndexTrans=zeros(length(traceYindexSmoothed),1);
xIndexTrans=yIndexTrans;
for n=1:length(traceYindexSmoothed)
    % x'=x*cos-y*sin
    % y'=x*sin+y*xos
    xIndexTrans(n)=xIndex(n)*cosd(theta)-traceYindexSmoothedCentered(n)*sind(theta);
    yIndexTrans(n)=xIndex(n)*sind(theta)+traceYindexSmoothedCentered(n)*cosd(theta);
    
    %or this way
    %xIndexTrans(n)=xIndex(n);
    %yIndexTrans(n)=traceYindexSmoothedCentered(n)-slope*xIndex(n);
end

%%
%going to trim ends off of transformed data before I try to fit to sine
xTrim=xIndexTrans(floor(length(xIndex)/4):length(xIndex)-floor(length(xIndex)/4));
yTrim=yIndexTrans(floor(length(xIndex)/4):length(xIndex)-floor(length(xIndex)/4));

%find zero crossings - maybe do this before we trim ends off
indexZeroXing=zeros(length(yTrim),1);%saves indeces to the left of each zero Xing
cnt=1;
for n=1:length(yTrim)-1
    if (yTrim(n+1)/yTrim(n)<0)
        indexZeroXing(cnt)=n;
        cnt=cnt+1;
    end
end

indexZeroXing=indexZeroXing(1:cnt-1);
guessP=2*mean(diff(indexZeroXing));
guessA=(max(yTrim)-min(yTrim))/2;
guessPhi=indexZeroXing(1);
guessY0=0;

% customSinFit=fittype('A*sin((2*pi*x/P)-phi)+y0','coefficients',{'A','P','phi','y0'});
customSinFit=fittype('A*sin(2*pi/P*(x-phi))+y0','coefficients',{'A','P','phi','y0'});
sinModel=fit(xTrim,yTrim,customSinFit,'startpoint',[guessA,guessP,guessPhi,guessY0]);
%ySinFit=sinModel.A*sin(2*pi/sinModel.P*(xTrim-sinModel.phi))+sinModel.y0;

sinParams=[sinModel.A, sinModel.phi];

end

