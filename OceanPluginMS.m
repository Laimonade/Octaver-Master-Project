% Ocean Plugin Matlab Style: Implemented with AnalysisAndSynthesisBuffer
classdef OceanPluginMS < audioPlugin 
    
    % User-defined properties to be used in the GUI
    properties
        Octave = 'Off'
        dwMix = 0.5
    end

    properties
        
        % Property determined by GUI 'octave'. Sets octave in semitone
        pitchFactor = 1;

        FFTReal;
        FFTImag;
        scaledReal;
        scaledImag;
        scaledSig;

        sampleRate = 44100;
            
        synthesisFactor = 2; % m * fftSize 
        frameNum = 0;
        
        %FrameTime Analysis period (s)
        %   Period of analysis frame in seconds, specified as a positive
        %   scalar between 0.02 and 0.04.
        %   The default is 0.03 seconds. This property is tunable.
        fftSize = 1024;  % 2048 @ <= 48k, 4096 @ 96k, 8192 @ 192k

        fftSize2;

        FrameTime = 0.02322;  % fftSize / sampleRate
        
        %Overlap Analysis frame overlap (%)
        %   Analysis frame overlap, specified in the range 0% to 80%.
        %   The default is 50%. This property is tunable.
        Overlap  = 75;
        
        % Gives a hope factor e.g: 50% = 2, 75% = 4
        hopFactor = 4 %100 / (100 - Overlap);

        %WindowType Type of analysis window
        %   Type of analysis window, specified as 'Triangle', 'Hann',
        %   'Hamming', 'Blackman', or 'Rectangle'.
        %   The default is 'Rectangle'. This property is tunable.
        WindowType = audioexample.WindowEnum.Hann;
         
        % Set action type for the overlap part (average or add)
        SynthesisOverlapAction = 'add';
        
        % Set samplerate for the A&S Buffer, modify 'sampleRate' if needed.
        SampleRate = 44100; % sampleRate;

    end
    
    properties(Constant, Access=private)
        TWOPI = 6.2831853071795864;
    end
    
    properties (Access = private)
        Buffer
    end

    properties (Constant)
%         PluginInterface = audioPluginInterface(...
%             'PluginName','OceanPluginMS', ...
%             'BackgroundImage', audiopluginexample.private.mwatlogo);

        PluginInterface = audioPluginInterface(...
            audioPluginParameter('Octave', 'Mapping', {'enum', '-2', '-1', 'Off', '1', '2'}), ...
            audioPluginParameter('WindowType',...
                'DisplayName','Analysis Window',...
                'Mapping',{'enum','Rectangle','Triangle','Hann','Hamming','Blackman'}),...
            audioPluginParameter('FrameTime',...
                'DisplayName','Analysis Frame',...
                'Label','s',...
                'Mapping',{'lin',0.01,0.04}),...
            audioPluginParameter('Overlap',...
                'DisplayName','Frame Overlap',...
                'Label','%',...
                'Mapping',{'int',0,80}));
    end
    

    methods
        function plugin = OceanPluginMS            
%         function setupImpl(plugin, ~)
            max_frame_len = 8192;
            plugin.FFTReal = zeros(max_frame_len, 1);
            plugin.FFTImag = zeros(max_frame_len, 1);
            plugin.scaledReal = zeros(max_frame_len, 1);
            plugin.scaledImag = zeros(max_frame_len, 1);
            plugin.scaledSig = zeros(max_frame_len, 1);

            sampleInput = zeros(1,1);
            plugin.Buffer = audiopluginexample.private.AnalysisAndSynthesisBuffer(sampleInput,44100,...
                                    'FrameTime',0.02322,...
                                    'AnalysisFrameOverlap',50,...
                                    'SynthesisFrameOverlap',50,...
                                    'WindowType',  audioexample.WindowEnum.Hann, ...
                                    'SynthesisOverlapAction', 'add');
            
            % Get useful info for analysis/synthesis
            plugin.fftSize2 = plugin.fftSize / plugin.hopFactor;
        end
        

 
        function out = process(plugin, in)
                       
            disp(length(in))
            disp(size(in))

            % Sum to mono for harmonic analysis
            monoIn = double(in); % Ensure doubles for analysis
            monoIn = plugin.sumToMono(monoIn);
            
            if plugin.Octave == "Off"
                out = in;

            else
            
                % Write to input buffer
                writeAnalysisBuffer(plugin.Buffer, monoIn);
                
                while (plugin.Buffer.AnalysisNumHopsAvailable > 0)
    
                    analysisFrame = readAnalysisBuffer(plugin.Buffer);
                    
                    % Compute FFT of windowed frame
                    currentFFT = complex(fft(analysisFrame, plugin.fftSize));
                    
                    % Compute multiplier 
                    multiplier = -(plugin.TWOPI * single(plugin.frameNum) / single(plugin.hopFactor * plugin.synthesisFactor * plugin.fftSize));
                    
    %                 semitoneFactor = octaveToPitch(plugin);
    
                    for bin = 1:plugin.fftSize
                        % Computing new bin number
%                         newBin = ceil(single(plugin.synthFactor * bin) * plugin.pitchFactor + 0.5);

                                                % m * k * a
                        newBin = round(plugin.synthesisFactor * plugin.pitchFactor * bin);

                                                    % analysis FFT size * a
                        if newBin >= 0 && newBin <= (plugin.fftSize * plugin.synthesisFactor)
%                             rl = real(currentFFT(bin));
%                             im = imag(currentFFT(bin));
%                             oscVar = single(newBin - (plugin.synthesisFactor * bin)) * multiplier;
%                             plugin.scaledReal(newBin) = (rl * cos(oscVar)) - (im * sin(oscVar));
%                             plugin.scaledImag(newBin) = (rl * sin(oscVar)) + (im * cos(oscVar));
                              plugin.scaledSig(newBin) = currentFFT(bin) * exp(1i * (newBin - (plugin.synthesisFactor * bin) * multiplier));
                        end
    
                    end
                    
                    % Adding back scaled real and imaginary part
%                     scaledSpectrum = complex(plugin.scaledReal, plugin.scaledImag);
                    
                    scaledSpectrum = complex(plugin.scaledSig);
                    
                    scaledSpectrum = [scaledSpectrum; conj(flip(scaledSpectrum(2:end-1)))]; % Reconstruct negative freqs.

                    % Taking Inverse FFT
                    scaledSignal = real(ifft(scaledSpectrum, (4 * plugin.fftSize)));    
                    
                    % Overlap add the --> New larger array has incompatible
                    % size with too small of a window
    
                    writeSynthesisBuffer(plugin.Buffer, scaledSignal(1:1024));
                    
                    plugin.frameNum = plugin.frameNum + 1;
    
                end
                
                plugin.frameNum = 1;
    
    
                % Check available unread data in output buffer
                NumUnreadSamples = checkSynthesisBufferMemory(plugin.Buffer);
                
                % Read from output buffer
                if (NumUnreadSamples > size(in,1)) && ...
                        (NumUnreadSamples > plugin.Buffer.AnalysisSamplesPerFrame)
                    y = readSynthesisBuffer(plugin.Buffer,size(in,1));
                    out = y(:,1);
                else
                    out = zeros(size(in));
                end
            end

% 
%             if plugin.frameNum == (plugin.hopFactor * plugin.synthFactor) 
%                 plugin.frameNum = 0;
%             end

            %%%% Perform pitch shift %%%
%             out = scaledSignal;
        
        end

        function reset(plugin)
            %reset Reset internal states to initial conditions
            % reset(noiseReducer) returns all private properties to their
            % initial conditions. In a DAW environment, reset is called
            % between uses or if the environment sample rate changes. In
            % MATLAB, call reset before an audio stream loop to mimic the
            % DAW environment.
            %
            plugin.Buffer.SampleRate = getSampleRate(plugin);
            reset(plugin.Buffer);
        end


    end
    
    %----------------------------------------------------------------------
    % PUBLIC METHODS : Listen for parameter tuning and update dependent
    % properties
    %----------------------------------------------------------------------
    methods
        %%%% SETUP A&S BUFFER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.FrameTime(plugin,val)
            plugin.Buffer.FrameTime = val; %#ok<MCSUP>
        end

        function val = get.FrameTime(plugin)
            val = plugin.Buffer.FrameTime;
        end

        function set.Overlap(plugin,val)
            plugin.Buffer.AnalysisFrameOverlap = val; %#ok<MCSUP>
            plugin.Buffer.SynthesisFrameOverlap = val; %#ok<MCSUP>
        end

        function val = get.Overlap(plugin)
            val = plugin.Buffer.AnalysisFrameOverlap;
        end

        function set.SampleRate(plugin, val)
            plugin.SampleRate = val;
        end
            
        function set.WindowType(plugin,val)
            switch val
                case audioexample.WindowEnum.Triangle
                    plugin.Buffer.WindowType = audioexample.WindowEnum.Triangle; %#ok<MCSUP>
                case audioexample.WindowEnum.Rectangle
                    plugin.Buffer.WindowType = audioexample.WindowEnum.Rectangle; %#ok<MCSUP>
                case audioexample.WindowEnum.Hann
                    plugin.Buffer.WindowType = audioexample.WindowEnum.Hann; %#ok<MCSUP>
                case audioexample.WindowEnum.Hamming
                    plugin.Buffer.WindowType = audioexample.WindowEnum.Hamming; %#ok<MCSUP>
                case audioexample.WindowEnum.Blackman
                    plugin.Buffer.WindowType = audioexample.WindowEnum.Blackman; %#ok<MCSUP>
            end
        end

        function val = get.WindowType(plugin)
            val = plugin.Buffer.WindowType;
        end
        
        function set.SynthesisOverlapAction(plugin, val)
            plugin.SynthesisOverlapAction = val;
        end

        function val = get.SynthesisOverlapAction(plugin)
            val = plugin.Buffer.SynthesisOverlapAction;
        end
        
        function set.Octave(plugin, val)
            plugin.Octave = val;
            switch val
                case '-2'
                    % plugin.pitchFactor = -24; %#ok<MCSUP>
                    plugin.pitchFactor = 0.25; %#ok<MCSUP>
                case '-1'
                    % plugin.pitchFactor = -12; %#ok<MCSUP>
                    plugin.pitchFactor = 0.5; %#ok<MCSUP>
                case 'Off'
                    plugin.pitchFactor = 1; %#ok<MCSUP>
                case '1'
                    % plugin.pitchFactor = 12; %#ok<MCSUP>
                    plugin.pitchFactor = 2; %#ok<MCSUP>
                case '2'
                    % plugin.pitchFactor = 24; %#ok<MCSUP>
                    plugin.pitchFactor = 4; %#ok<MCSUP>
            end
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods(Access = private)

        function monoOut = sumToMono(~,in)
            [~,n] = size(in);
            
            if n == 2
                monoOut = sum(in, 2);
            end
        end
    end
end
