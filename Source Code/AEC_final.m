%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Digital Signal Processing                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          As a Part of Digital Signal Processing Laboratory              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Acoustic Echo Cancellation with Double Talk Detection            %
%       Done By: 1) Suraj Khan (22455660)                                 %
%                2) Shrishti Saha Shetu (22464279)                        %
%                3) Raza Azam ()                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%set filter length
N = 2001;

[IR]=load('IR.Mat', 'IR');
IR = (IR.IR);


%set step size
ALPHA = 1;

mem_time = 0.125;   % seconds of memory for correlator
freeze_time = 3.5;    % Min of 2 Seconds of filter turn off for pauses in between speech
del_processing_time= 0.125; % Filtering delayed

%load signals
[x, fs] = audioread('ref.wav');
[d, fs] = audioread('mic_sig_spe_20.wav');
IR= IR*fs;

%Have a persistence for the the impulse response over time
len=length(x);
h_str=zeros(len,2001);
distance=zeros(len,1);
%run NLMS & DD
[e, h_hat, d_hat,h_str] = NLMS_DD(x, d, N, ALPHA, fs, mem_time, freeze_time, del_processing_time,h_str);
%System Performance Measure

%System distance
for i=1:len
IR_e=IR - h_str(i,:).';
distance(i)=((IR_e')*IR_e)/(IR'*IR);
end

x_axis=(1:length(x))/16000;



%SDM
figure
plot(x_axis,10*log10(distance));
title('System Distance Measure')
ylabel('System Distance (Logairthmic Scale)');
xlabel('Time in Seconds');

%Microphone and Error Signal
figure
subplot(211); plot(x_axis,d);title('Microphone Signal'); xlabel('Time in Seconds'); ylabel('Amplitude');
subplot(212); plot(x_axis,e);title('Error Signal'); xlabel('Time in Seconds'); ylabel('Amplitude');

%ERLE
ERLE = 10 * log10(filter(0.1, [1 -0.98], d.^2) ./ filter(0.1, [1 -0.98], e.^2));
figure
plot(x_axis(1:end-4000),ERLE(1:end-4000)); xlabel('Time in Seconds'); ylabel('ERLE [dB]');
title('ERLE Plot');

%Impulse response
figure;
subplot(211); plot(IR);
title('True Impulse Response'); 
xlabel('Filter Taps'); 
ylabel('Coefficient Value');
subplot(212); 
plot(h_hat);
title('Estimated Impulse Response'); 
xlabel('Filter Taps'); 
ylabel('Coefficient Value');

%listen to the residual
player_test=audioplayer(e, fs);