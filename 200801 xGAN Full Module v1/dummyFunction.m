function [time,voltage] = dummyFunction(t_res, f_rep, waveform, rawWaveformTime, rawWaveformSignal)
%DUMMYFUNCTION This function selects one of several output waveforms to
%plot then outputs the function and fourier transform of that function
%   Functions include
%       - 'custom' - A custom input
%       - 'gPulse' - A Gaussian Pulse
%       - 'gEdge' - A Gaussian Edge
%       - 'RC_Edge' - An RC Edge
%       - 'Cos' - A cosine function
%   Internal parameters are set by directly modifying the function

%% Sets axes
if(strcmp(waveform,'custom')
    time = rawWaveformTime(1:end-1)  %Removing last time point
    n = size(time)
    n = n(2)
    t_res  = time(2)-time(1); % temporal resolution
    f_res = 1/t_res;
    t_rep = interpTimeDbl(end) -time(1); %Rep rate shifted up to real value
    f_rep = 1/t_rep;
    f_LT_array = linspace(-f_res/2,f_res/2,n)
else
    % t_res = 0.05*c.nano; % temporal resolution
    f_res = 1/t_res;
    % f_rep = 10*c.mega; %10MHz repetition rate
    t_rep = 1/f_rep; %period of the repetition

    timeN = -t_rep/2:t_res:-t_res;          % From negative start just to before 0
    time0P = 0:t_res:t_rep/2-t_res;         % From 0 to maximum just before maximum positive
    time = [timeN time0P];                  % Combine two

    tLength = length(time);                  %Should be divisible by 2!! for halves of response
    n = tLength;


    fP = f_res*(0:(tLength/2))/tLength;         % Since FT symmetric now
    %f  = f_res*(-tLength/2:tLength/2)/tLength;        % Full FT frequency domain
    f = linspace(-f_res/2, f_res/2, n)
end
%% Defines functions (Modify Coefficients)
switch waveform
    case 'gEdge'
        gEdgeCoef = 10e15;
        gEdgeN = exp(-gEdgeCoef*timeN.^2);
        gEdge0P = (gEdgeN(1)+1-exp(-gEdgeCoef*timeN.^2));
        fn = [gEdgeN gEdge0P];
    case 'gPulse'
        gaussCoef = 5e15;
        fn = exp(-gaussCoef*time.^2);
        
    case 'RC_Edge'
        rcCoef = 1e-9;
        RC_EdgeN = 1-exp(-((timeN-timeN(1))/rcCoef));
        RC_Edge0P  = RC_EdgeN(end)*exp(-time0P/rcCoef);
        fn = [RC_EdgeN RC_Edge0P];
        fn = rescale(fn,-0.5, 0.5);
        
    case 'Cos'
        nPer = 50; % number of periods
        fn = cos(2*nPer*pi*f_rep*time);

    otherwise
        disp('invalid input')
end

end
