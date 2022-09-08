% A class performing the sliding DFT in real-time. Computes a DFT and its
% inverse on a sample-by-sample basis instead of the traditional
% block-by-block approach. This allows to obtain an output every sample and
% in some configuration, reducing the latency to 1 sample. 

classdef OPTI_RT_SDFT < handle
    properties (Access = private)

        % Circular buffer containing input samples
        inBuffer = [];
        
        % SDFT size (Assigned through object init)
        dftSize {mustBeNumeric}

        % Read pointer keeping track of where we are in input buffer
        % Wrapping around buffer every N samples
        readPointer = 1;

        % A & S twiddle factors
        twiddleAnalysis = complex([]);
        twiddleSynthesis = complex([]);

        % Array storing frequency domain values pre frequency domain
        % windowing. 
        freqsArray = complex([]);
        
        % Damping factor r and it's r^N power. Helps stability 
        % https://www.intechopen.com/chapters/54042
        r;
        rPowN;

        % Weight to compensate windowing in frequency domain
        windowWeight;

        % Weight to correct scaling of inverse SDFT
        isdftWeight;

        % Flag keeping track of data readiness 
        % (Init latency of N samples while all the init 0's are replaced)
        % NOTE: This is only an initial latency while the DFT starts up and
        % does not reflect latency once running. 
        dataLengthFlag = false;
    end

    % Init private constants for the plugin
    properties(Constant, Access=private)
        TWOPI = 6.2831853071795864;
    end

    properties (Access = public)
        % Array storing frequency domain values post frequency domain
        % windowing. 
        dft = [];

        % Damping factor to help reduce common stability issues due to loss
        % of precision of the twiddle factors over time. Issues due to pole
        % displacement outside the unit circle
        dampingFactor = 1; % --> to check
        timer;
        idx;
    end

    methods

        % Constructor
        function obj = OPTI_RT_SDFT(SDFT_SIZE)

            % Make sure DFT length is passed to constructor
            if nargin == 1
                obj.dftSize = SDFT_SIZE;
            else
                error("You need to specify the DFT length as argument");
            end
            
            % Using proper DFT notation for length
            N = obj.dftSize;
            
            % Resize & init buffers
            obj.inBuffer = zeros(N, 1);
            obj.freqsArray = complex(zeros(N, 1));
            obj.dft = complex(zeros(N, 1));
           
            % Resize twiddles
            obj.twiddleAnalysis = complex(zeros(N,1));
            obj.twiddleSynthesis = complex(zeros(N,1));

             % Update the DFT
            obj.r = obj.dampingFactor;
            obj.rPowN = power(obj.r, N);

            % Weight to compensate windowing in frequency domain
            obj.windowWeight = 1/obj.dftSize;

            obj.isdftWeight = 1 / (N * 0.5);

            % Compute analysis twiddle factor
            for n = 0:N-1
                factor = obj.TWOPI * n / N;
                obj.twiddleAnalysis(n+1) = exp(1i * factor);
            end
            
            % Compute synthesis twiddle factor
            for t = 0:N-1
                % We have the choice to choose from which value we
                % reconstruct the signal in TD. If we start with the first
                % value of the DFT, we save on computation however we have
                % 1 window size of delay. Taking the last value of the DFT
                % has 0 latency but is computationally heavy as we have to
                % perform N computations. Usual value is N/2 offering a
                % good comprise between latency and computation load.

                START_VALUE = N;
                factor = -obj.TWOPI * t * START_VALUE / N;
                obj.twiddleSynthesis(t+1) = exp(1i * factor);
            end
           
            obj.timer = zeros(1000,1);
            obj.idx =1;
        end




        % Update the calculation with a new sample
        % Returns true if the data are valid (because enough samples have been
        % presented), or false if the data are invalid.
        function validFlag = sdft(obj, new_sample)

            % Save current sample for this loop's calculation
            previous_sample = obj.inBuffer(obj.readPointer);
            
            % Store new incoming sample
            obj.inBuffer(obj.readPointer) = new_sample;

            % Computing DFT shifting all current DFT points and adding the
            % new sample
            obj.freqsArray = obj.twiddleAnalysis(1:obj.dftSize) .* (obj.r .* obj.freqsArray - obj.rPowN * previous_sample + new_sample);
            
            % Apply the Hanning window using its frequency response and multiplying it with the signal- ~0.018sec
            % Acts as a filter after the DFT step
            obj.dft = simd_windowing_mex(obj.freqsArray, obj.dft, obj.dftSize, obj.windowWeight);

            % Increment read pointer to grab next sample in input array
            obj.readPointer = obj.readPointer + 1;
            
            % Check that pointer has not fallen off the edge of vector and
            % not trying to access random memory
            if obj.readPointer > obj.dftSize
                obj.dataLengthFlag = true;
                obj.readPointer = 1;
            end
            
            % After 1 loop through (N samples), we set the data flag to
            % true
            validFlag = obj.dataLengthFlag;
        end
    

        function tdSample = isdft(obj, dft)
            
            tdAccumulator = 0;
            
            % ISDFT_TYPE: last sample (1 sample latency but expensive computation)
            tdAccumulator = tdAccumulator + sum(real(dft .* obj.twiddleSynthesis));           
            
            tdSample = tdAccumulator * obj.isdftWeight;

        end
    

    end
end