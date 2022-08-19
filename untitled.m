% overlap * zeroPad;
% roots = getRoots(4,4);
% cmpx = complex(-20.6, 7.65);  
% disp(roots(2))
% disp(cmpx)
% 
% disp(cmpx *roots(2))

evalFreq = 1000.0; % Hz
sampleRate = 44100.0; % Hz
blockSize = 1; % FFT size
zeropad = 1; % Used to increase synthesis DFT reoslution
pitchShiftRatio = 2; % Double in frequency

error = expectedError(evalFreq, sampleRate, blockSize, zeropad, pitchShiftRatio);
disp(["Expected error: " error])