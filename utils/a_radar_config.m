% Configuration refer to TI mmWave Studio
% TI Param -- Aways check before data process
numADCSample = 1.280000e+02; 
adcSampleRate = 2.280000e+06; %Hz/s 
startFreqConst = 7.700000e+10; %Hz 
chirpSlope = 6.334300e+13; %Hz/s 
chirpIdleTime = 9.700000e-05; %s 
adcStartTimeConst = 6.000000e-06; %s 
chirpRampEndTime = 6.314000e-05; %s 
NumDevices = 4; 
framePeriodicty = 1.000000e-01; 
numChirpsInLoop = 1.200000e+01; %s 
nchirp_loops = 50; 
numTxAnt = 12; 

%% constants
speedOfLight = 3e8;
Fc = 77e9; % GHz
c = 3e8;

%% derived parameters
DopplerFFTSize = 2^(ceil(log2(nchirp_loops)));
numChirpsPerFrame = nchirp_loops*numChirpsInLoop;%nchirp_loops*numTxAnt;
chirpRampTime       = numADCSample/adcSampleRate;
chirpBandwidth      = chirpSlope * chirpRampTime; % Hz
chirpInterval       = chirpRampEndTime + chirpIdleTime;
carrierFrequency    = startFreqConst +  (adcStartTimeConst + chirpRampTime/2)*chirpSlope; % Hz center frequency
lambda              = speedOfLight/carrierFrequency;
maximumVelocity     = lambda / (chirpInterval*4) ; % m/s
maxRange            = speedOfLight*adcSampleRate*chirpRampTime/(2*chirpBandwidth);
numSamplePerChirp   = round(chirpRampTime*adcSampleRate);
rangeFFTSize        = 2^(ceil(log2(numSamplePerChirp)));
numChirpsPerVirAnt  = nchirp_loops;
rangeResolution     = speedOfLight/2/chirpBandwidth;
rangeBinSize        = rangeResolution*numSamplePerChirp/rangeFFTSize;
velocityResolution  = lambda / (2*nchirp_loops * chirpInterval*numTxAnt);
velocityBinSize     = velocityResolution*numChirpsPerVirAnt/DopplerFFTSize;


RANGE_FFT = numADCSample;
AZIM_FFT = 128; % fft across antennas -- equal to the BP_fft_test angle_NFFT
MAX_RANGE = RANGE_FFT * rangeResolution; % msize()sizeter RANGE_FFT * range_bin_resolution

second_fft = imgs;
[num_of_frames, num_of_range_bins, num_of_angle_bins] = size(second_fft);

% correcting fs
fs = 1/framePeriodicty*nchirp_loops;

% Creating axis labels
maxdist = MAX_RANGE;
range_NFFT = RANGE_FFT;
angle_NFFT = AZIM_FFT;
doppler_NFFT = 64;
sine_theta = -2*linspace(-0.5,0.5,angle_NFFT);
cos_theta = sqrt(1-sine_theta.^2);
[range, sine_theta_mat] = ndgrid(linspace(0,maxdist,range_NFFT),sine_theta);
[~, cos_theta_mat] = ndgrid(linspace(0,maxdist,range_NFFT),cos_theta);
x_axis = range.*cos_theta_mat;
y_axis = range.*sine_theta_mat;

