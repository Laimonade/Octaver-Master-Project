function [outputVector] = ocean(inData, winSize, hopFactor, scaleFactor, synthFactor)

%% Parameters
fftFrameSize2 = fftFrameSize/2;
hopSize = fftFrameSize/hopFactor;
inFifoLatency = fftFrameSize-stepSize;
outFrameFullSize = fftFrameSize * synthFactor;

% Pitch scaling factor
alpha = 2^(scaleFactor/12);

% Intermediate constants
hopOut = round(alpha * hop);

% Hanning window for overlap-add
wn = hann(winSize*2+1);
wn = wn(2:2:end);

%% Source file

x = inData(:,1);

% Rotate if needed
if size(x,1) < size(x,2)
   
    x = transpose(x);

end

% Left zero padding
x = [zeros(hop*3,1) ; x];

%% Initialization

% Create a frame matrix for the current input
% Analysis overlap add
[y,numberFramesInput] = createFrames(x,hop,winSize);

% Create a frame matrix to receive processed frames
numberFramesOutput = numberFramesInput;
outputy = zeros(numberFramesOutput,winSize);

% Initialize cumulative phase
phaseCumulative = 0;

% Initialize previous frame phase
previousPhase = 0;

for index=1:numberFramesInput

%% Analysis
    
    % Get current frame to be processed
    currentFrame = y(index,:);
    
    % Window the frame
    currentFrameWindowed = currentFrame .* wn' / sqrt(((winSize/hopSize)/2));
    
    % Get the FFT: returns an array of complex values (real, img) (mag,
    % phase)
    currentFrameWindowedFFT = fft(currentFrameWindowed);
    
    % Get the magnitude (3)
    magFrame = abs(currentFrameWindowedFFT);
    
    % Get the angle (phase, the phase is wrapped) (4)
    phaseFrame = angle(currentFrameWindowedFFT);
    
%% Processing    
        
    % Get the phase difference between 2 consecutive frames
    % This phase difference lies between -pi and pi (5)
    deltaPhi = phaseFrame - previousPhase;
    previousPhase = phaseFrame;
    
    % Remove the expected phase difference (5)
    deltaPhiPrime = deltaPhi - hop * 2*pi*(0:(winSize-1))/winSize;
    
    % Map to -pi/pi range, could use wrapToPi(deltaPhiPrimeMod) (5)
    deltaPhiPrimeMod = mod(deltaPhiPrime+pi, 2*pi) - pi;

    % Get the true frequency (6)
    trueFreq = 2*pi*(0:(winSize-1))/winSize + deltaPhiPrimeMod/hop;
    
    % Get the final phase (7) (8)
    phaseCumulative = phaseCumulative + hopOut * trueFreq;    
    
    % Remove the 60 Hz noise. This is not done for now but could be
    % achieved by setting some bins to zero.
   
%% Synthesis    
    
    % Get the magnitude
    outputMag = magFrame;
    
    % Produce output frame
    outputFrame = real(ifft(outputMag .* exp(1i*phaseCumulative)));
     
    % Save frame that has been processed
    outputy(index,:) = outputFrame .* wn' / sqrt(((winSize/hopOut)/2));
        
end

%% Finalize

% Overlap add in a vector
outputTimeStretched = fusionFrames(outputy,hopOut);

% Resample with linearinterpolation
outputTime = interp1((0:(length(outputTimeStretched)-1)),outputTimeStretched,(0:alpha:(length(outputTimeStretched)-1)),'linear');

% Return the result
outputVector = outputTime;

return