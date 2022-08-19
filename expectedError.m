function error = expectedError (evalFreq, sampleRate, blockSize, zeropad, pitchShiftRatio)
    fa = evalFreq;
    fs = sampleRate;
    N = blockSize;
    m = zeropad;
    k = pitchShiftRatio;

    error = ( ((fa + (fs / (2 * m * N))) * pitchShiftRatio) + (fs / (2 * m * N)) ) / (k * fa);
end