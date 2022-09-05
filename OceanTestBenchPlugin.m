% This code is a test bench audio plugin implementing the Juillerat Ocean algorithm for pitch
% shifting, allowing to tweak various parameters. 
% This code is inspired by the implementation that can be found
% on https://www.pitchtech.ch/ and the paper https://ieeexplore.ieee.org/document/5685027?arnumber=5685027.
% This plugin performs pitch shifting by computing the STFT of an audio stream, repositioning the frequency
% content to a new bin whose location is dependant on the pitch factor, and
% finally adjusting the phase of each bin to maintain vertical phase
% coherence. This pitch shifting system is supported by a Weighted Overlap
% Add (WOLA) infrastructure. The demodulation part is a work in progress as
% it does not seem to work properly. 

classdef OceanTestBenchPlugin < audioPlugin & matlab.System

    % User-defined properties to be used in the GUI
    % Init properties that the end-user interacts with.
    properties
        Octave = 'Off'
        Zeropad = '1'
        FFTLen = '4096'
        Demodulation = 'Off'
        Overlap = '4'
        DryWetMix = 1.0
    end

    % Init properties that the end-user does not interact with directly.
    properties (Access = private)

        % Property determined by GUI 'octave'. Sets octave in semitone
        pitchShiftRatio;

        blockSize; % 'N' in the paper
        overlap; % 'O' in the paper
        hopSize; % size of input and output audio chunks, N/O
        zeroPad; % 'm' in the paper
        nbInputBins; % Number of spectral bins, input
        nbOutputBins; % Number of spectral bins, output

        processHopIndex; % Running index of processed STFT frames since the beginning
        outputHopIndex; % Running index of synthesized STFT frames since the beginning

        % Buffers
        windowedInputBlock;
        windowedOutputBlock;
        analysisWindow;
        synthesisWindow;
        inputSpectrum;
        unityRoots;
        outputSpectrum;
        cycleLength;
        scaledRePart;
        scaledImPart;
        outputComplex;
        overflowWarned;

        % Size of input & ouput circular buffer
        bufferSize;

        % Circular input buffer
        inputBuffer;
        inputBufferPointer;
        hopCounter;

        % Circular output buffer
        outputBuffer;
        outputBufferInitLatency;
        outputBufferWritePointer;
        outputBufferReadPointer;

    end

    % Init private constants for the plugin
    properties(Constant, Access=private)
        TWOPI = 6.2831853071795864;
    end


    properties (Constant)

        PluginInterface = audioPluginInterface(...
            audioPluginParameter('Octave', 'Mapping', {'enum', '-2', '-1', 'Off', '1', '2'}), ...
            audioPluginParameter('FFTLen', 'Mapping', {'enum', '512', '1024', '2048', '4096'}), ...
            audioPluginParameter('Overlap', 'Mapping', {'enum', '1', '2', '4', '8'}), ...
            audioPluginParameter('Zeropad', 'Mapping', {'enum', '1', '2', '4'}), ...
            audioPluginParameter('Demodulation', 'Mapping', {'enum', 'On', 'Off'}), ...
            audioPluginParameter('DryWetMix', 'Mapping',{'lin',0,1}));
    end


    methods (Access = protected)

        function setupImpl(plugin, ~)

            %%%%%% INIT USEFUL VARIABLES %%%%%%
            plugin.blockSize = str2double(plugin.FFTLen); % 'N' in the paper
            plugin.overlap = str2double(plugin.Overlap); % 'O' in the paper
            plugin.hopSize = plugin.blockSize / plugin.overlap; % size of input and output audio chunks, N/O
            plugin.zeroPad = 1; % 'm' in the paper
            plugin.nbInputBins = (plugin.blockSize / 2) + 1; % Number of spectral bins, input
            plugin.nbOutputBins = ((plugin.blockSize * plugin.zeroPad) / 2) + 1; % Number of spectral bins, output
            plugin.processHopIndex = -(plugin.overlap - 1); % Running index of processed STFT frames since the beginning
            plugin.outputHopIndex = plugin.processHopIndex -1; % Running index of synthesized STFT frames since the beginning


            %%%%%% INIT BUFFERS %%%%%%

            plugin.bufferSize = 16384;

            % Initialize input buffer
            plugin.inputBuffer = zeros(plugin.bufferSize, 1, 'single');
            plugin.inputBufferPointer = 1;
            plugin.hopCounter = 1;

            % Initialize output buffer
            plugin.outputBuffer = zeros(plugin.bufferSize * plugin.zeroPad, 1, 'single');
            plugin.outputBufferInitLatency = plugin.hopSize ;
            plugin.outputBufferWritePointer = plugin.outputBufferInitLatency;
            plugin.outputBufferReadPointer = 1;


            % Init windowed input array
            plugin.windowedInputBlock = zeros(plugin.blockSize, 1, 'single');
            plugin.windowedOutputBlock = zeros(plugin.blockSize * plugin.zeroPad, 1, 'single');


            % Init analysis & synthesis windows
            plugin.analysisWindow = single(hann(plugin.blockSize));
            plugin.synthesisWindow = plugin.analysisWindow;

            % Init complex input & output spectrum
            plugin.inputSpectrum = complex(zeros(plugin.nbInputBins, 1), 0);
            plugin.outputSpectrum = complex(zeros(plugin.nbOutputBins, 1), 0);

            % Init scaled real and imaginary part
            plugin.scaledRePart = zeros(plugin.nbOutputBins, 1, 'single');
            plugin.scaledImPart = zeros(plugin.nbOutputBins, 1, 'single');

            plugin.cycleLength = plugin.overlap * plugin.zeroPad;

            plugin.unityRoots = complex(getRoots(plugin.cycleLength, plugin.cycleLength));
            plugin.overflowWarned = false;

        end


        % Main process loop
        function out = stepImpl(plugin, in)
            if length(in) ~= plugin.hopSize
                err = ['Input buffer needs to be the same size as hop size: ', num2str(plugin.hopSize), ' samples.'];
                error(err)
            end

            % Convert to float and sum to mono
            monoIn = plugin.sumToMono(single(in));
            
            % TO UNCOMMENT IF WE WANT TRUE BYPASS
            % If octaver not used, bypass expensive computation
            %             if plugin.Octave == "Off"
            %                 out = monoIn;
            %             else


            % Wrap around buffer once end is reached
            if plugin.inputBufferPointer + length(monoIn) - 1 > plugin.bufferSize
                plugin.inputBufferPointer = 1;
            end

            % Add system in buffer to circular inputBuffer
            plugin.inputBuffer(plugin.inputBufferPointer : plugin.inputBufferPointer ...
                + length(monoIn) - 1) = monoIn;

            % Increment circular inputBuffer pointer by one hop
            plugin.inputBufferPointer = plugin.inputBufferPointer + length(monoIn);

            % Get the OLA output samples
            output = plugin.outputBuffer(plugin.outputBufferReadPointer : plugin.outputBufferReadPointer ...
                + length(monoIn) - 1, 1);

            % Clear last output samples
            plugin.outputBuffer(plugin.outputBufferReadPointer : plugin.outputBufferReadPointer ...
                + length(monoIn) - 1, 1) = zeros(length(monoIn), 1, 'single');

            % Scale output down by the overlap factor
            output = output .* single(1 / (0.75* plugin.overlap));

            % Increment circular inputBuffer pointer
            plugin.outputBufferReadPointer = int32(plugin.outputBufferReadPointer + length(monoIn));

            % Increment output read pointer
            % Wrap around buffer once end is reached
            if plugin.outputBufferReadPointer + length(monoIn) - 1 > plugin.bufferSize
                plugin.outputBufferReadPointer = 1;
            end

            plugin.hopCounter = plugin.hopCounter + length(monoIn);

            % Because we force audioIn = hopSize, this condition is
            % always true.
            if plugin.hopCounter >= plugin.hopSize
                plugin.hopCounter = 1;

                % Do FFT
                unwrappedBuffer = zeros(plugin.blockSize, 1, 'single');

                % Using modulo to find the index in matlab leads to an
                % issue (which took me ages to figure out). The entire
                % point to use the modulo operator is to know when we
                % have to wrap around the buffer. If we have a buffer
                % of 1024 samples, when we get to the 1024 sample, we
                % have 1024 % 1024 = 0. This works well for languages
                % that start their index at 0 but that's not the case
                % in matlab. We are also missing the last element of
                % the buffer as in matlab we need [1:1024] and not
                % [0:1023]. This leads to 'add 1' to the length of the
                % buffer size. Because we are increasing the buffer by
                % 256 values, we are wrapping around bufferSize/blockSize,
                % leading to sometimes wrap around between [769-1024:1-256].
                % When this happen, we add 1 to the index to scale the
                % i to go from [0:255] to [1:256] --> 'indexTrick'.

                indexTrick = 0;
                for i = 1:plugin.blockSize

                    index = mod((plugin.inputBufferPointer + i - plugin.blockSize + plugin.bufferSize), ...
                        plugin.bufferSize + 1);

                    if index ~= 0
                        inCircBuffIndex = index + indexTrick;
                    else
                        inCircBuffIndex = 1;
                        indexTrick = 1;
                    end

                    unwrappedBuffer(i) = plugin.inputBuffer(inCircBuffIndex);
                end

                % Apply analysis window using vectorisation (element-wise mult)
                plugin.windowedInputBlock = unwrappedBuffer .* plugin.analysisWindow;

                % Perform FFT on the windowed input block
                plugin.inputSpectrum = fft(plugin.windowedInputBlock, plugin.blockSize);

                % Shift signal by the specified octave -
                % modifies plugin.outputSpectrum in place

                % If octaver not used, bypass expensive computation
                if plugin.Octave == "Off"
                    reconstrutedAndShiftedSignal = plugin.inputSpectrum;
                else
                    plugin.octaveShiftSignal();


                    % Rebuild the conjugate part of the signal not used during
                    % the FFT.
                    reconstrutedAndShiftedSignal = complex(vertcat(plugin.outputSpectrum, conj(flip(plugin.outputSpectrum(2:end-1)))));
                end
                % Perform IFFT of pitch shifted signal
                modSignal = real(ifft(reconstrutedAndShiftedSignal, plugin.blockSize * plugin.zeroPad));

                % Crop signal to blocl size and apply window output block with synthesis window
                plugin.windowedOutputBlock = modSignal(1:plugin.blockSize) .* plugin.synthesisWindow;

                indexTrick = 0;
                for i = 1:plugin.blockSize

                    index = mod((plugin.outputBufferWritePointer + i - 1), plugin.bufferSize + 1);

                    if index ~= 0
                        outCircBuffIndex = index + indexTrick;
                    else
                        outCircBuffIndex = 1;
                        indexTrick = 1;
                    end

                    if plugin.Demodulation == "Off"
                        demodulatorFactor = 1;
                    else
                        demodulatorFactor = plugin.demodulator(i, plugin.outputHopIndex);

                    end

                    plugin.outputBuffer(outCircBuffIndex, 1) = (plugin.outputBuffer(outCircBuffIndex) + plugin.windowedOutputBlock(i))...
                        * demodulatorFactor;
                end

                % Updated output buffer write pointer to start at the
                % next hop
                outputPointer =  mod((plugin.outputBufferWritePointer + plugin.hopSize), plugin.bufferSize);
                if outputPointer ~= 0
                    plugin.outputBufferWritePointer = mod((plugin.outputBufferWritePointer + plugin.hopSize), plugin.bufferSize);
                else
                    plugin.outputBufferWritePointer = 1;
                end

            end

            out = (output * plugin.DryWetMix) + ((1 - plugin.DryWetMix) * monoIn);

            plugin.processHopIndex = plugin.processHopIndex + 1;
            plugin.outputHopIndex = plugin.outputHopIndex + 1;

        end



        function resetImpl(plugin)
            % reset Reset internal states to initial conditions
            % reset(OceanOlaPlugin) returns all private properties to their
            % initial conditions. In a DAW environment, reset is called
            % between uses or if the environment sample rate changes. In
            % MATLAB, call reset before an audio stream loop to mimic the
            % DAW environment.

%             reset(OceanOlaPlugin);
            % reset input buffer
            plugin.inputBuffer = zeros(plugin.bufferSize, 1, 'single');
            plugin.inputBufferPointer = 1;
            plugin.hopCounter = 1;

            % reset output buffer
            plugin.outputBuffer = zeros(plugin.bufferSize * plugin.zeroPad, 1, 'single');
            plugin.outputBufferWritePointer = plugin.outputBufferInitLatency;
            plugin.outputBufferReadPointer = 1;

            % reset windowed input array
            plugin.windowedInputBlock = zeros(plugin.blockSize, 1, 'single');
            plugin.windowedOutputBlock = zeros(plugin.blockSize * plugin.zeroPad, 1, 'single');

            % reset analysis & synthesis windows
            plugin.analysisWindow = single(hann(plugin.blockSize));
            plugin.synthesisWindow = plugin.analysisWindow;

            % reset complex input & output spectrum
            plugin.inputSpectrum = complex(zeros(plugin.nbInputBins, 1), 0);
            plugin.outputSpectrum = complex(zeros(plugin.nbOutputBins, 1), 0);

            % reset scaled real and imaginary part
            plugin.scaledRePart = zeros(plugin.nbOutputBins, 1, 'single');
            plugin.scaledImPart = zeros(plugin.nbOutputBins, 1, 'single');

            plugin.cycleLength = plugin.overlap * plugin.zeroPad;

            plugin.hopSize = plugin.blockSize / plugin.overlap;

        end
    end



    methods

        function monoOut = sumToMono(~,in)
            [~,n] = size(in);

            if n == 2
                monoOut = sum(in, 2);
            end
        end

        function octaveShiftSignal(plugin)

            inputFFT = plugin.inputSpectrum;

            % Copy DC
            plugin.outputSpectrum(1) = inputFFT(1);

            % Clear previous bins
            plugin.outputSpectrum(2:end) = zeros(plugin.nbOutputBins - 1, 1);

            % TODO: Replace cycle by frame (p) to respect paper convention
            plugin.cycleLength = plugin.overlap * plugin.zeroPad;
            cycleIndex = mod((plugin.processHopIndex + plugin.cycleLength * 2), plugin.cycleLength);

            multiplier = single(-plugin.TWOPI * cycleIndex / plugin.cycleLength);
            
             
            for bin = 2:plugin.nbInputBins
                % Compute actual bin taking account of padding
                paddedBin = bin * plugin.zeroPad;

                % Compute new shifted bin
                newBin = floor((paddedBin * plugin.pitchShiftRatio) + 0.5); % b:= floor(m*k*a + 0.5)
               
                % Make sure bin is in between range of outputBin
                if newBin > 1 && newBin <= plugin.nbOutputBins 
                    realPart = real(inputFFT(bin));
                    imPart = imag(inputFFT(bin));

                    % full theta equation
                    phaseMult = single(newBin - (paddedBin) * multiplier);

                    scaledRePrt = (realPart * cos(phaseMult)) - (imPart * sin(phaseMult));
                    scaledImPrt = (realPart * sin(phaseMult)) + (imPart * cos(phaseMult));
                    plugin.outputComplex = complex(scaledRePrt, scaledImPrt);

                    plugin.outputSpectrum(newBin) = plugin.outputComplex;

                end
            end
        end
        

        % Finds the exact modulation window which the pitch shifting
        % process creates. The aim is to demodulate the output by dividing
        % the it by this modulation window. 
        function demodulationWindow = demodulator(plugin, index, hopIndex)
            result = single(0.0);

            for i = 1:plugin.overlap
                % Maybe try i - 1 & index - 1
                offset = mod(((i - 1) * plugin.hopSize) + index, plugin.blockSize);
                if offset == 0
                    offset = 1;
                end

                pitchShiftratioForFrame = plugin.pitchShiftRatio;

                result = result + plugin.getModifiedAnalysisWindow(offset, hopIndex, pitchShiftratioForFrame)...
                    * plugin.synthesisWindow(offset);
            end

            threshold = single(0.1);

            if result <= threshold

                if ~plugin.overflowWarned
                    disp("Cannot accurately demodulate. Please increase overlap!");
                    plugin.overflowWarned = true;

                end

                demodulationWindow =  1.0 / threshold;

            else
                demodulationWindow =  1.0 / result;

            end

        end

        function modWindow = getModifiedAnalysisWindow(plugin, offsetIndex, hopIndex, pitchShiftRatio)

            % hopIndex % cycleLength (Euclidian modulo, assuming hopIndex >= -cycleLength*2)
            cycleIndex = int32(mod((hopIndex + plugin.cycleLength * 2), plugin.cycleLength));

            % Calculate 'pitchShiftRatio - zeroPad', modulo 'cycleLength', Euclidian modulo
            % With no zero padding (zeroPad == 1), this is 'pitchShiftRatio - 1', modulo 'overlap'
            psrMinusPad = mod((pitchShiftRatio + (plugin.cycleLength - plugin.zeroPad)), plugin.cycleLength);

            shift = int32(( plugin.blockSize * cycleIndex * psrMinusPad) / plugin.cycleLength); % s * N

            offset = mod((offsetIndex * pitchShiftRatio / plugin.zeroPad + shift), plugin.blockSize); % (k*t+s*N) % N

            if offset ~= 0
                modWindow = plugin.analysisWindow(offset);
            else
                modWindow = plugin.analysisWindow(1);
            end

        end


    end

    %----------------------------------------------------------------------
    % PUBLIC METHODS : Listen for parameter tuning and update dependent
    % properties
    %----------------------------------------------------------------------
    methods
        %%%% Getters & Setters methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% E.G Of A Getter Method %%%
        % function val = get.SynthesisOverlapAction(plugin)
        %   val = plugin.Buffer.SynthesisOverlapAction;
        % end

        function set.Octave(plugin, val)
            plugin.Octave = val;
            switch val
                case '-2'
                    plugin.pitchShiftRatio = 0.25; %#ok<MCSUP>
                case '-1'
                    plugin.pitchShiftRatio = 0.5; %#ok<MCSUP>
                case 'Off'
                    plugin.pitchShiftRatio = 1; %#ok<MCSUP>
                case '1'
                    plugin.pitchShiftRatio = 2; %#ok<MCSUP>
                case '2'
                    plugin.pitchShiftRatio = 4; %#ok<MCSUP>
            end

            resetImpl(OceanTestBenchPlugin)

        end

        function set.Zeropad(plugin, val)
            plugin.Zeropad = val;
            plugin.zeroPad = str2double(val); %#ok<MCSUP>
            plugin.resetImpl;
        end

        function set.Demodulation(plugin, val)
            plugin.Demodulation = val;
            plugin.resetImpl;
        end

        function set.DryWetMix(plugin, val)
            plugin.DryWetMix = val;
        end

        function set.Overlap(plugin, val)
            plugin.Overlap = val;
            plugin.overlap = str2double(val); %#ok<MCSUP>
            plugin.resetImpl;
        end

        function set.FFTLen(plugin, val)
            plugin.FFTLen = val;
            plugin.blockSize = str2double(val); %#ok<MCSUP>
            plugin.resetImpl;
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

