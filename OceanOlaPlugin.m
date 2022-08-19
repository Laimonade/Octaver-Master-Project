% OverlapAddTest
classdef OceanOlaPlugin < audioPlugin & matlab.System
    
    % User-defined properties to be used in the GUI
    % Init properties that the end-user interacts with.
    properties
        Octave = 'Off'
        Zeropad = '1'
        dwMix = 0.5
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
            audioPluginParameter('Zeropad', 'Mapping', {'enum', '1', '2', '4'}));
    end
    

    methods (Access = protected)

        function setupImpl(plugin, ~)

            %%%%%% INIT USEFUL VARIABLES %%%%%%
            plugin.blockSize = 1024; % 'N' in the paper
            plugin.overlap = 4; % 'O' in the paper
            plugin.hopSize = plugin.blockSize / plugin.overlap; % size of input and output audio chunks, N/O
%             plugin.zeroPad = 2; % 'm' in the paper
            plugin.nbInputBins = (plugin.blockSize / 2) + 1; % Number of spectral bins, input
            plugin.nbOutputBins = ((plugin.blockSize * plugin.zeroPad) / 2) + 1; % Number of spectral bins, output
            plugin.processHopIndex = -(plugin.overlap - 1); % Running index of processed STFT frames since the beginning
            plugin.outputHopIndex = plugin.processHopIndex; % Running index of synthesized STFT frames since the beginning


            %%%%%% INIT BUFFERS %%%%%%
            
            plugin.bufferSize = int32(16384);

            % Initialize input buffer
            % plugin.inputBlock = dsp.AsyncBuffer;
            % write(plugin.inputBlock, zeros(plugin.blockSize, 1, 'single'));
            plugin.inputBuffer = zeros(plugin.bufferSize, 1, 'single');
            plugin.inputBufferPointer = 1;
            plugin.hopCounter = 1;
            
            % Initialize output buffer
            % plugin.outputBlock = dsp.AsyncBuffer;
            % write(plugin.outputBlock, zeros(plugin.blockSize * plugin.zeroPad));
            plugin.outputBuffer = zeros(plugin.bufferSize * plugin.zeroPad, 1, 'single');
            plugin.outputBufferInitLatency = plugin.blockSize - plugin.hopSize; % plugin.hopSize;
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

        end
        

        % Main process loop
        function out = stepImpl(plugin, in)    
            if length(in) ~= plugin.hopSize
                error('Buffer needs to be the same size as hop size')
            end
            
            % Convert to float and sum to mono 
            monoIn = plugin.sumToMono(single(in));
            
            % If octaver not used, bypass expensive computation
            if plugin.Octave == "Off"
                out = monoIn;
            else
                
                % Wrap around buffer once end is reached
                if plugin.inputBufferPointer + length(monoIn) - 1 > plugin.bufferSize
                    plugin.inputBufferPointer = 1;    
                end

                % Add system in buffer to circular inputBuffer
                plugin.inputBuffer(plugin.inputBufferPointer : plugin.inputBufferPointer ...
                    + length(monoIn) - 1) = monoIn;

                % Increment circular inputBuffer pointer
                plugin.inputBufferPointer = plugin.inputBufferPointer + length(monoIn);
                
                % Get the OLA output samples 
                output = plugin.outputBuffer(plugin.outputBufferReadPointer : plugin.outputBufferReadPointer ...
                    + length(monoIn) - 1, 1);
                
                % Clear last output samples
                plugin.outputBuffer(plugin.outputBufferReadPointer : plugin.outputBufferReadPointer ...
                    + length(monoIn) - 1, 1) = zeros(length(monoIn), 1, 'single');
                
                % Scale output down by the overlap factor
                % output = output .* single(1 / plugin.overlap);
                
                % Increment circular inputBuffer pointer
                plugin.outputBufferReadPointer = int32(plugin.outputBufferReadPointer + length(monoIn));
                
                % Increment output read pointer
                % Wrap around buffer once end is reached
                if plugin.outputBufferReadPointer + length(monoIn) - 1 > plugin.bufferSize
                    plugin.outputBufferReadPointer = 1;   
                end
                
                plugin.hopCounter = plugin.hopCounter + length(monoIn);

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
                    plugin.octaveShiftSignal();
                    
                    % Rebuild the conjugate part of the signal not used during
                    % the FFT.
                    reconstrutedAndShiftedSignal = complex(vertcat(plugin.outputSpectrum, conj(flip(plugin.outputSpectrum(2:end-1)))));
                    
                    % Perform IFFT of pitch shifted signal
                    modSignal = real(ifft(reconstrutedAndShiftedSignal, plugin.blockSize * plugin.zeroPad));
                    
                    % Window output block with synthesis window
                    plugin.windowedOutputBlock = modSignal(1:plugin.blockSize) .* plugin.synthesisWindow;
                    
                    indexTrick = 0;
                    for i = 1:plugin.blockSize
                        index = mod((plugin.outputBufferWritePointer + i), plugin.bufferSize + 1);
                        if index ~= 0
                            outCircBuffIndex = index + indexTrick;
                        else
                            outCircBuffIndex = 1;
                            indexTrick = 1;
                        end
                        plugin.outputBuffer(outCircBuffIndex, 1) = plugin.outputBuffer(outCircBuffIndex) + plugin.windowedOutputBlock(i);
                    end

                    % Updated output buffer write pointer to start at the
                    % next hop
                    plugin.outputBufferWritePointer = mod((plugin.outputBufferWritePointer + plugin.hopSize), plugin.bufferSize);
                end
                
                out = output;

            end

        end
        
        

        function resetImpl(plugin)
            %reset Reset internal states to initial conditions
            % reset(noiseReducer) returns all private properties to their
            % initial conditions. In a DAW environment, reset is called
            % between uses or if the environment sample rate changes. In
            % MATLAB, call reset before an audio stream loop to mimic the
            % DAW environment.
            
            reset(OceanOlaPlugin);

            plugin.inputBuffer = zeros(plugin.bufferSize, 1, 'single');
            plugin.inputBufferPointer = 1;
            plugin.hopCounter = 1;
            
            % Initialize output buffer
            % plugin.outputBlock = dsp.AsyncBuffer;
            % write(plugin.outputBlock, zeros(plugin.blockSize * plugin.zeroPad));
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
        end

        function set.Zeropad(plugin, val)
            plugin.Zeropad = val;
            plugin.zeroPad = str2double(val); %#ok<MCSUP>
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
            
            multiplier = single(-plugin.TWOPI / plugin.cycleLength);
             
            for bin = 2:plugin.nbInputBins
                % Compute actual bin taking account of padding
                paddedBin = bin * plugin.zeroPad;

                % Compute new shifted bin
                newBin = floor((paddedBin * plugin.pitchShiftRatio) + 0.5); % b:= floor(m*k*a + 0.5)
                
%                 PHASEMETHOD = "github";
                PHASEMETHOD = "paper";


                if PHASEMETHOD == "github"
                    % Make sure bin is in between range of outputBin
                    if newBin > 0 && newBin < plugin.nbOutputBins
                        rePart = real(inputFFT(bin));
                        imPart = imag(inputFFT(bin));
    
                        % Compute (b-ma) % mO from paper
                        if newBin >= paddedBin
                            cycleShift = mod((newBin - paddedBin), plugin.cycleLength);
                        else
                            cycleShift = mod(plugin.cycleLength - (paddedBin - newBin), plugin.cycleLength);
                        end
                        
                        % (b-ma)p % mO
                        phaseShift = mod((cycleShift * cycleIndex), plugin.cycleLength);
                        
                        % full theta equation
                        phaseMult = single(phaseShift * multiplier);
                        
                        if phaseShift ~= 0
                            scaledRePrt = (rePart * cos(phaseMult)) - (imPart * sin(phaseMult));
                            scaledImPrt = (imPart * sin(phaseMult)) + (imPart * cos(phaseMult));
                            plugin.outputComplex = complex(scaledRePrt, scaledImPrt);
                        end
                        plugin.outputSpectrum(newBin) = plugin.outputComplex;

                    end
                elseif PHASEMETHOD == "paper"
                   if newBin > 0 && newBin < plugin.nbOutputBins
                       work = complex(inputFFT(bin));

                       % Compute (b-ma) % mO from paper
                        if newBin >= paddedBin
                            cycleShift = mod((newBin - paddedBin), plugin.cycleLength);
                        else
                            cycleShift = mod(plugin.cycleLength - (paddedBin - newBin), plugin.cycleLength);
                        end
                        
                        % (b-ma)p % mO
                        phaseShift = mod((cycleShift * cycleIndex), plugin.cycleLength);
                        
                        if phaseShift ~= 0
                            index = mod((plugin.cycleLength - phaseShift), plugin.cycleLength);
                            if index ~= 0
                                plugin.outputComplex = work * plugin.unityRoots(index);
                            else
                                plugin.outputComplex = work * plugin.unityRoots(1);
                            end
                        end
                        plugin.outputSpectrum(newBin) = complex(plugin.outputComplex);

                   end

                end
               
                
            end
        end


    end
end

