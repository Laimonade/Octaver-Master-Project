% OverlapAddPluginTemplate
classdef OverlapAddPluginTemplate < audioPlugin & matlab.System
    
    % User-defined properties to be used in the GUI
    % Init properties that the end-user interacts with.
    properties
        
    end
    
    % Init properties that the end-user does not interact with directly.
    properties (Access = private)
        
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
        
        % Size of input & ouput circular buffer
        bufferSize;

        % Circular input buffer
        inputBuffer;
        inputBufferPointer;
        hopCounter;

        % Circular output buffer
        outputBuffer;
        outputBufferWritePointer;
        outputBufferReadPointer;
            
        
        
    end
    
    % Init private constants for the plugin
    properties(Constant, Access=private)

    end
    

    properties (Constant)

        PluginInterface = audioPluginInterface();
    end
    

    methods (Access = protected)

        function setupImpl(plugin, ~)

            %%%%%% INIT USEFUL VARIABLES %%%%%%
            plugin.blockSize = 1024; % 'N' in the paper
            plugin.overlap = 4; % 'O' in the paper
            plugin.hopSize = plugin.blockSize / plugin.overlap; % size of input and output audio chunks, N/O
            plugin.zeroPad = 1; % 'm' in the paper
            plugin.nbInputBins = (plugin.blockSize / 2) + 1; % Number of spectral bins, input
            plugin.nbOutputBins = ((plugin.blockSize * plugin.zeroPad) / 2) + 1; % Number of spectral bins, output
            plugin.processHopIndex = -(plugin.overlap - 1); % Running index of processed STFT frames since the beginning
            plugin.outputHopIndex = plugin.processHopIndex; % Running index of synthesized STFT frames since the beginning


            %%%%%% INIT BUFFERS %%%%%%
            
            plugin.bufferSize = int32(16384);
%             plugin.bufferSize = int32(1024);


            % Initialize input buffer
            % plugin.inputBlock = dsp.AsyncBuffer;
            % write(plugin.inputBlock, zeros(plugin.blockSize, 1, 'single'));
            plugin.inputBuffer = zeros(plugin.bufferSize, 1, 'single');
            plugin.inputBufferPointer = 1;
            plugin.hopCounter = 1;
            
            % Initialize output buffer
            % plugin.outputBlock = dsp.AsyncBuffer;
            % write(plugin.outputBlock, zeros(plugin.blockSize * plugin.zeroPad));
            plugin.outputBuffer = zeros(plugin.bufferSize, 1, 'single');
            plugin.outputBufferWritePointer = plugin.hopSize;
            plugin.outputBufferReadPointer = 1;

            % Initialize overlap add buffer
            % plugin.overlapAddBuffer = dsp.AsyncBuffer;
            % write(plugin.overlapAddBuffer, zeros(plugin.blockSize - plugin.hopSize, 1, 'single'));
            plugin.overlapAddBuffer = zeros(plugin.blockSize - plugin.hopSize, 1, 'single');
            
            
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
            
        end
        

        % Main process loop
        function out = stepImpl(plugin, in)    
%             if length(in) ~= plugin.hopSize
%                 error('Buffer needs to be the same size as hop size')
%             end
            
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
%                 output = output .* single(1 / plugin.overlap);
                
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
                    % buffer size. Because we are increasing the buffer by 256 values, 
                    % we are wrapping around bufferSize/blockSize, leading
                    % to sometimes wrap around between [769-1024:1-256].
                    % When this happen, we add 1 to the index to scale the
                    % i to go from [0:255] to [1:256] --> 'indexTrick'.

                    indexTrick = 0;
                    for i = 1:plugin.blockSize
                        md = mod((plugin.inputBufferPointer + i - plugin.blockSize + plugin.bufferSize), ...
                            plugin.bufferSize + 1); 

                        if md ~= 0
                            inCircBuffIndex = md + indexTrick;
                        else
                            inCircBuffIndex = 1;
                            indexTrick = 1;
                        end
                        unwrappedBuffer(i) = plugin.inputBuffer(inCircBuffIndex);
                    end
                    
                    windowedUnwrappedBuffer = unwrappedBuffer .* plugin.analysisWindow;

                    spectrum = fft(windowedUnwrappedBuffer, plugin.blockSize);
            
                    modSignal = ifft(spectrum, plugin.blockSize);
                    
                    windowedModSignal = modSignal .* plugin.synthesisWindow;
                    
                    indexTrick = 0;
                    for i = 1:plugin.blockSize
                        md = mod((plugin.outputBufferWritePointer + i), plugin.bufferSize + 1);
                        
                        if md ~= 0
                            outCircBuffIndex = md + indexTrick;
                        else
                            outCircBuffIndex = 1;
                            indexTrick = 1;
                        end
                        plugin.outputBuffer(outCircBuffIndex, 1) = plugin.outputBuffer(outCircBuffIndex) + windowedModSignal(i);
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
            plugin.num = 1;

            reset(OverlapAddPluginTemplate);

            plugin.inputBuffer = zeros(plugin.bufferSize, 1, 'single');
            plugin.inputBufferPointer = 1;
            plugin.hopCounter = 1;
            
            % Initialize output buffer
            % plugin.outputBlock = dsp.AsyncBuffer;
            % write(plugin.outputBlock, zeros(plugin.blockSize * plugin.zeroPad));
            plugin.outputBuffer = zeros(plugin.bufferSize, 1, 'single');
            plugin.outputBufferWritePointer = plugin.blockSize - plugin.hopSize;
            plugin.outputBufferReadPointer = 1;

            % reset overlap add buffer
            plugin.overlapAddBuffer = zeros(plugin.blockSize - plugin.hopSize, 1, 'single');
           
            
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
        
        %%%% E>G Of A Setter Method
%         function set.Octave(plugin, val)
%             switch val
%                 case '-2'
%                     plugin.pitchShiftRatio = 0.25; %#ok<MCSUP>
%                 ...
              

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    methods

        function monoOut = sumToMono(~,in)
            [~,n] = size(in);
            
            if n == 2
                monoOut = sum(in, 2);
            end
        end
        


    end
end
