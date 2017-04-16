% Script for generating sounds filtered using the HRTFs
% THIS SCRIPT IS ALSO USED FOR 
% get the HRIR file
% Azimuth conversion: for us, 0 is straight ahead
% LEFT  is n (negative) 1-90, 1-90 in HRIR coordinates
% RIGHT is p (positive) 1-90, 359-270 in HRIR coordinates

clear
clc
%Enter list of subjects(IRCAM database)
subjectID = 1021;  %change this for different settings  %1023, 1021, 1005, 1008

hrirFileDir = ['/Users/lab/Data/HRIR(Ircam)/IRC_' num2str(subjectID)...
                    '/RAW/MAT/HRIR/'];
load([hrirFileDir '/IRC_' num2str(subjectID) '_R_HRIR.mat']);


hrir_pts  = size(l_hrir_S.content_m,2); % number of sample in HRIR
az = l_hrir_S.azim_v(l_hrir_S.elev_v==0);
naz  = length(az);
sample_rate_hz = 44100;

%************************* FILTERED NOISE *****************
% some standard params for generating sounds
pulse_reps = 5;
pulse_gap = 0.01;
sound_len_sec = 0.05;   
ramp_len_sec   = 0.005;   %sounds start and end with cosine ramps
ramp_samples   = round(ramp_len_sec * sample_rate_hz);
bpfilter       = design(fdesign.bandpass('N,F3DB1,F3DB2',...
                                 10,100,15e3,sample_rate_hz));

%Sounds are band-pass filtered Gaussian waveforms with ramps
%at the beginning and the end to avoid "splatter"
sound_samples = round(sound_len_sec*sample_rate_hz);
sound_base    = filter(bpfilter, randn(sound_samples, 1));
sound_base(1:ramp_samples) = sound_base(1:ramp_samples,:).*...
    cos(linspace(3*pi/2,2*pi,ramp_samples))';
sound_base(end-ramp_samples+1:end) = sound_base(end-ramp_samples+1:end).*...
    cos(linspace(0,pi/2,ramp_samples))';

hrir_az  = rem(360-az, 360); % remap degrees for hrir indexing
ssz  = (sound_samples + round(0.1*sound_samples))*pulse_reps;
sdat = nan(ssz, 2, naz);
aux = zeros(round(pulse_gap*sample_rate_hz),1);
sound_base = [sound_base ; aux];
sound_base = repmat(sound_base, pulse_reps,1);
%**********************************************************




% %*********************** Frequency sweep *****
% hrir_az  = rem(360-az, 360);
% sound_len_sec = 0.2;   
% sound_samples = round(sound_len_sec*sample_rate_hz);
% t = linspace(0,sound_len_sec,sound_samples);
% w1 = 2*pi*100;
% w2 = 2*pi*10000;
% a1 = w1*sound_len_sec/(log(w2/w1));
% sound_base = sin(a1*(exp(log(w2/w1).*t./sound_len_sec)-1));
% 
% %band pass filter the sound
% bpfilter       = design(fdesign.bandpass('N,F3DB1,F3DB2',...
%                                  10,100,15e3,sample_rate_hz));
% ramp_len_sec   = 0.005;   %sounds start and end with cosine ramps
% ramp_samples   = round(ramp_len_sec * sample_rate_hz);
% sound_base    = filter(bpfilter, sound_base);
% 
% ssz  = sound_samples + round(0.1*sound_samples);
% sdat = nan(ssz, 2, naz);
% sound_base  = sound_base';
% %*********************************************


for aa = 1:naz
   index = find(l_hrir_S.elev_v==0 & l_hrir_S.azim_v==hrir_az(aa), 1);
   cLeft = conv(sound_base, l_hrir_S.content_m(index,:)');
   cRight = conv(sound_base, r_hrir_S.content_m(index,:)');
   len = min(ssz,length(cLeft));
   cS = [cLeft(1:len) cRight(1:len)];
   sdat(1:len,:,aa) = (cS/max(abs(cS(:))))/1.01;
end


zeroIdx = find(az==0);
ahead = audioplayer(sdat(:,:,zeroIdx), sample_rate_hz);


for aa = 1:naz  
    disp(az(aa))
    aps = audioplayer(sdat(:,:,aa), sample_rate_hz);
    %sound(sdat(:,:,aa), sample_rate_hz);
    pause(0.4);
    %playblocking(ahead);
    playblocking(aps);
end





