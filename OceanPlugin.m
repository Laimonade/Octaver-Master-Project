classdef OceanPlugin < audioPlugin 
    
    properties
        pitchFactor = 1
        fftSize = 1024
        windowSize = 256
        dwMix = 0.5
    end
    
    properties(Constant, Access=private)
        TWOPI = 6.2831853071795864;
    end
    
    properties(Access = private, Hidden)
        inputBuffer;
        outputBuffer;
        outputAccumulatorBuffer;
        analysisBuffer;
        FFTReal;
        FFTImag;
        scaledReal;
        scaledImag;

        nFFT = 2048;
        hopFactor = 2;
        synthFactor = 4;
        frameNum = 0;
        wn;

    end

    properties (Constant)
        pluginInterface = audioPluginInterface(...
            audioPluginParameter('pitchFactor', 'Mapping', {'lin', -12, 12}))
    end
    

    methods
        function plugin = OceanPlugin
             % Initialize buffer for harmonic analysis
            plugin.analysisBuffer = dsp.AsyncBuffer;
            write(plugin.analysisBuffer, [0; 0]);
            read(plugin.analysisBuffer, 2);
            
            plugin.inputBuffer = dsp.AsyncBuffer;
            write(plugin.inputBuffer, [0 0; 0 0]);
            read(plugin.inputBuffer);
            
            plugin.outputBuffer = dsp.AsyncBuffer;
            write(plugin.outputBuffer, [0 0; 0 0]);
            read(plugin.outputBuffer);
            
            plugin.outputAccumulatorBuffer = dsp.AsyncBuffer;
            write(plugin.outputAccumulatorBuffer, [0 0; 0 0]);
            read(plugin.outputAccumulatorBuffer);
              

            % Hanning window for overlap-add
            win = hann(plugin.nFFT*2+1);
            plugin.wn = win(2:2:end);

            max_frame_len = 8192;
            FFTReal = zeros(max_frame_len, 1);
            FFTImag = zeros(max_frame_len, 1);
            scaledReal = zeros(max_frame_len, 1);
            scaledImag = zeros(max_frame_len, 1);
            write(plugin.outputAccumulatorBuffer, zeros(max_frame_len * 2, 2));
            write(plugin.analysisBuffer, zeros(plugin.nFFT));


    end
 
        function out = process(plugin, in)
            
            
            % Get useful info for analysis/synthesis
            fs = getSampleRate(plugin);
            n_fft = plugin.nFFT; % 2048 @ <= 48k, 4096 @ 96k, 8192 @ 192k
            n_fft2 = n_fft / plugin.hopFactor;

            % Sum to mono for harmonic analysis
            monoIn = double(in); % Ensure doubles for analysis
            monoIn = plugin.sumToMono(monoIn);
            
            % Write to buffer
            write(plugin.analysisBuffer, monoIn);
            
            % If analysis buffer is full enough (here fft), perform fft
            if plugin.analysisBuffer.NumUnreadSamples >= n_fft

                % Async buffer, num of rows, overlap size
                analysisFrame = read(plugin.analysisBuffer, n_fft, n_fft2);
                
                % Apply hann windowing
                analysisWinFrame = analysisFrame .* plugin.wn' / sqrt(((plugin.nFFT/plugin.hopFactor)/2));
                
                % Compute FFT of windowed frame
                currentFFT = fft(analysisWinFrame);
                
                % Compute multiplier 
                multiplier = -(plugin.TWOPI * single(plugin.frameNum) / single(plugin.hopFactor * plugin.synthFactor * n_fft));

                for bin = 1:n_fft2
                    % Computing new bin number
                    newBin = ceil(single(plugin.synthFactor * bin) * plugin.pitchFactor + 0.5);

                    if newBin >= 0 && newBin <= (n_fft * plugin.synthFactor)
                        rl = real(currentFFT(bin));
                        im = imag(currentFFT(bin));
                        oscVar = single(newBin - (plugin.synthFactor * bin)) * multiplier;
                        plugin.scaledReal(newBin) = (rl * cos(oscVar)) - (im * sin(oscVar));
                        plugin.scaledImag(newBin) = (rl * sin(oscVar)) + (im * cos(oscVar));
                    end

                end
                
                % Adding back scaled real and imaginary part
                scaledSpectrum = plugin.scaledReal + plugin.scaledImag;
                
                % Taking Inverse FFT
                scaledSignal = ifft(scaledSpectrum);
                
                % Overlap add the --> New larger array has incompatible
                % size with too small of a window
                synthesisWinFrame = scaledSignal .* plugin.wn' / sqrt(((winSize/hop)/2));

          

            end

            plugin.frameNum = plugin.frameNum + 1;

            if plugin.frameNum == (plugin.hopFactor * plugin.synthFactor) 
                plugin.frameNum = 0;
            end

            %%%% Perform pitch shift %%%
            out = synthesisWinFrame;
        
        end


    end

    methods(Access = private)

        function monoOut = sumToMono(~,in)
            [~,n] = size(in);
            
            if n == 2
                monoOut = sum(in, 2);
            end
        end
    end
end
