%%
%% add zff and data path
path_code = 'C:\Research_drive_RSP\whispp\Whispp_Works-20220704T074841Z-001\Whispp_Works\';
addpath C:\Research_drive_RSP\whispp\Whispp_Works-20220704T074841Z-001\Features_From_ZFR_and_ZTW\Features_From_ZFR_and_ZTW\
%addpath C:\Research_drive_RSP\whispp\Whispp_Works-20220704T074841Z-001\Whispp_Works\cleanPathologicalSpeech\
%addpath C:\Research_drive_RSP\whispp\Whispp_Works-20220704T074841Z-001\Whispp_Works\cleanWhisperSpeechWithBabble\
addpath( genpath( 'C:\Research_drive_RSP\whispp\Whispp_Works-20220704T074841Z-001\Whispp_Works\modifying_zeroFreqFiltering\' ) ) ;
path_clean = 'C:\Research_drive_RSP\whispp\Whispp_Works-20220704T074841Z-001\pathological_speech_no_backround_noise\'; 
path_path = 'C:\Research_drive_RSP\whispp\Whispp_Works-20220704T074841Z-001\whispered_speech_with_babble\';
addpath( path_clean , path_path ); 
%%
%%
%% load the audio 
filename = [path_path '20220524_Binnentuin_restaurant_single_channel_iPhone_whispp_mic.wav'];
[x , f] = audioread(filename); 
if(f~=16000); x = resample(x,16000,f); f = 16000; end
lx = length(x); 
xs = x;
% xs = x(6:11.5*f); 
%%% lowpass and highpass filter
disp('performing bandpass')
[b_lpf,a_lpf] = butter( 4 , 0.1 , 'low' );
xlf = filter( b_lpf , a_lpf , xs );
[b_hpf,a_hpf] = butter( 4 , 0.01  , 'high' ); % 0.05
xlf = filter( b_hpf , a_hpf , xlf );
xlf = xlf/max(abs(xlf)); 

[b_hpf,a_hpf] = butter( 4 , 0.2  , 'high' ); % 0.05
xhpf = filter( b_hpf , a_hpf , xs );
xhpf = xhpf/max(abs(xhpf)); 

y = xlf+xhpf; y = y/max(abs(y)); 
% sound( [xs ; y ] ,f );

% strN = strsplit(filename,'\'); 
% the_file = ['BPF\' char(strN(end))] ; 
% audiowrite(the_file, y ,f); 
% disp( ' done writing audio ' )
% clear strN xhpf xlpf 
%% %% -- %% PLOTTING FIGURES
[ss_xs]=createSpectrogram( xs , 20*f/1000 , 10*f/1000 , 1024 ); %ss_xs = averageSpgram( ss_xs , 16 );
[ss_xhf]=createSpectrogram( xhpf , 20*f/1000 , 10*f/1000 , 1024 ); %ss_xlf = averageSpgram( ss_xlf , 16 );
[ss_xlf]=createSpectrogram( xlf , 20*f/1000 , 10*f/1000 , 1024 ); %ss_xlf = averageSpgram( ss_xlf , 16 );
[ss_y]=createSpectrogram( y , 20*f/1000 , 10*f/1000 , 1024 ); %ss_xlf = averageSpgram( ss_xlf , 16 );

% ss_xhf = averageSpgram( ss_xhf , 7 );
% ss_xlf = averageSpgram( ss_xlf , 7 );
% ss_y = averageSpgram( ss_y , 15 );
fax = f/1024:f/1024:f/2; 
[m,n] = size(ss_xs); 

clf;
subplot(211); 
surf( 1:n , fax , ss_xlf ); colormap( flipud( gray ) ); shading interp ; view([0 90]); ylim([0 3000]); %view([-90 80]);

subplot(212); 
surf( 1:n , fax , ss_xhf ); colormap( flipud( gray ) ); shading interp ; view([0 90]); ylim([0 3000]); %view([-90 80]);

%view([-90 80]);
%% NMF based analysis
[ss_xhf]=createSpectrogram( xhpf , 20*f/1000 , 10*f/1000 , 1024 ); %ss_xlf = averageSpgram( ss_xlf , 16 );
[ss_xlf]=createSpectrogram( xlf , 20*f/1000 , 10*f/1000 , 1024 ); %ss_xlf = averageSpgram( ss_xlf , 16 );
ss_xlf = min(min(ss_xlf))/10 + removeNAN(ss_xlf); 
ss_xhf = min(min(ss_xhf))/10 + removeNAN(ss_xhf); 

[w1,h1] = nnmf( ss_xlf , 10 ); 
[w2,h2] = nnmf( ss_xhf , 10 ); 

subplot(211); plot(w1); grid on; 
subplot(212); plot(w2); grid on; 

%% modified zero frequency filtering 


%% linear prediction 
windSize = 30 * f/1000 ;
windShift = 10 * f/1000 ;
lporder = 16 ;
preemphasis = 1 ; 
sig = xs( 1e3 : end ) ; 
sig = sig /max( abs ( sig ) ) ; 
fs = f ; 
[lpres,lpcoeff,sig]=lpAnalysis(sig,fs,windSize,windShift,lporder,preemphasis);
lpresA = [ diff(lpres) ; 0 ];

lpresB = zfsig( sig , fs , 2 );
lpresC = zfsig( sig , fs , .5 );

lpres = lpresA/max(abs(lpresA)) + lpresB/max(abs(lpresB)) + lpresC/max(abs(lpresC)) ; 
% lpres = lpresA + lpresB + lpresC ;
lpres = lpres./max(abs(lpres)); 
[syn_sig]=lpSynthesis(lpcoeff,lpres,fs,windSize,windShift) ;

% mat_B = createSpectrogram( lpresB , 30*f/1000 , 10*f/1000 , 1024 ) ;
% mat_D = createSpectrogram( lpresD , 30*f/1000 , 10*f/1000 , 1024 ) ;
% 
% clf;
% subplot(211); surf( mat_B ); colormap( flipud( gray ) ); shading interp ; view([0 90]); ylim([0 112]); xlim([100 600]);
% subplot(212); surf( mat_D ); colormap( flipud( gray ) ); shading interp ; view([0 90]); ylim([0 112]); xlim([100 600]);

% [lpc_mat]=createSpectrogram_lpc( xs , 40*f/1000 , 10*f/1000 , 1024 , order );
% A = log(lpc_mat); A = lpc_mat; 
% A = averageSpgram( lpc_mat , 6 );
% figure; surf( A ); colormap( flipud( gray ) ); shading interp ; view([0 90]); ylim([0 212]); xlim([200 600]);


%% spectral average and reconstruction
winL = 20*f/1000 ; 
winS = 10*f/1000 ;


bx = buffer( xs , winL , winL - winS , 'nodelay' );
[m,n] = size(bx); 
dft_mat = zeros( winL , n ); 
for ix = 1 : n
    fx = fft( bx(:,ix) );
    dft_mat(:,ix) =  fx; 
end
T = angle(dft_mat);
ss_av = averageSpgram( abs(dft_mat) , 5 );
ss_av = ss_av.*( cos(T)+1i*sin(T) );


yy_xs = overlapNadd( ss_av , winL , winS ); 
sound(real(yy_xs),f)
yy = real(yy_xs); yy = yy / max(yy) ;
strN = strsplit(filename,'\'); 
% plot( real(yy_xs) ); grid; %ylim([-1.5 1.5]);
the_file = ['avSpec\' char(strN(end))];
audiowrite(the_file, yy ,f); 
disp( ' done writing audio 2' )
%%
%%

%%