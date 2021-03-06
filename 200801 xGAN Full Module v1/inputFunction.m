function [time, signal, freq, fftSignal] = inputFunction(waveform, t_res, f_rep,  rawWaveformTime, rawWaveformSignal, tuningParam, timeShift)
%DUMMYFUNCTION This function selects one of several output waveforms to
%plot then outputs the function and fourier transform of that function
%   Functions include
%       - 'Custom' - A custom input
%       - 'gPulse' - A Gaussian Pulse
%       - 'gEdge' - A Gaussian Edge
%       - 'RC_Edge' - An RC Edge
%       - 'Cos' - A cosine function
%   Internal parameters are set by directly modifying the function

%% Sets time/freq axes
if(strcmp(waveform,'Custom'))
    time = rawWaveformTime(1:end);  %Removing last time point
    n = size(time);
    n = n(2);
    t_res  = time(2)-time(1); % temporal resolution
    f_res = 1/t_res;
    t_rep = rawWaveformTime(end) -time(1); %Rep rate shifted up to real value
    f_rep = 1/t_rep;
    freq = linspace(-f_res/2,f_res/2,n);
    
else
    % t_res = 0.05*c.nano; % temporal resolution
    f_res = 1/t_res;
    % f_rep = 10*c.mega; %10MHz repetition rate
    t_rep = 1/f_rep; %period of the repetition

    timeN = -t_rep/2:t_res:-t_res;          % From negative start just to before 0
    time0P = 0:t_res:t_rep/2;         % From 0 to maximum just before maximum positive
    time = [timeN time0P];                  % Combine two

    tLength = length(time);                  %Should be divisible by 2!! for halves of response
    n = tLength;
    
    %fP = f_res*(0:(tLength/2))/tLength;         % Since FT symmetric now
    %f  = f_res*(-tLength/2:tLength/2)/tLength;        % Full FT frequency domain
    freq = linspace(-f_res/2, f_res/2, n);
end
%% Defines functions (Modify Coefficients)
switch waveform
    case 'gEdge'
        gEdgeCoef = tuningParam;
        gEdgeN = -exp(-gEdgeCoef*timeN.^2);
        gEdge0P = -(gEdgeN(1)+1-exp(-gEdgeCoef*flip(time0P).^2));
        fn = flip([gEdgeN gEdge0P]);
        
    case 'gPulse'
        gaussCoef = tuningParam;
        fn = exp(-gaussCoef*time.^2);
        
    case 'RC_Edge'
        rcCoef = tuningParam;
        RC_EdgeN = 1-exp(-((timeN-timeN(1))/rcCoef));
        RC_Edge0P  = RC_EdgeN(end)*exp(-time0P/rcCoef);
        fn = [RC_EdgeN RC_Edge0P];
        fn = rescale(fn,-0.5, 0.5);
        
    case 'Heaviside'
        delay = 50;
        lengthOff = tuningParam;
        ltot = n;
        heav.A = ones(1,delay);
        heav.B = zeros(1,lengthOff);
        heav.C = ones(1,ltot - delay - lengthOff);

        fn = [heav.A heav.B heav.C];
        
    case 'ErrorFn'
        slowEdge = 3;
        fastLength = 200;
        testFnPos = erf(linspace(-tuningParam, tuningParam, fastLength));
        testFnNeg = erf(linspace(-slowEdge, slowEdge, n - fastLength));
        testFnPos = flip(testFnPos);
        fn = [testFnNeg testFnPos];
    
    case 'GaussFilt'
        fftFilter = exp(-(tuningParam*linspace(-1, 1, 2001).^2));
        fftFilter = fft(fftFilter);
        
        A = 200;
        B = 100;
        lin = 800;
        
        heav.A = ones(1,A); 
        heav.B = zeros(1,B);
        heav.lin = linspace(0, 1, lin);
        heav.C = ones(1,2001-(A+B+lin));
        fnIn = [heav.A heav.B heav.lin heav.C];
        fnMix = fft(fnIn).*fftFilter;
        fn = ifft(fnMix);
        
        
    case 'Cos'
        nPer = tuningParam; % number of periods
        fn = cos(2*nPer*pi*f_rep*time);

    case 'Custom'
        fn = rawWaveformSignal(1:end);
    otherwise
        disp('invalid input')
end
signal = circshift(fn, timeShift);

fftSignal = fft(signal); 
end

