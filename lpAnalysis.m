function[lpres,lpcoeff,sig]=lpAnalysis(sig,fs,windSize,windShift,lporder,preemphasis)

%
%function[lpres,lpcoeff,sig8k]=lpAnalysis(sig,fs,windSize,windShift,lporder,preemphasis)
%
%
% Performs LP analysis of the speech signal
%
%

%
% Coded by: Anand Joseph, SVL, IIIT Hyderabad
% Revision: 1a, 15/07/2009
%

% if fs>8000
%     sig8k=resample(sig,8000,fs);
%     sig=sig8k;
% end

windSize = round(windSize*fs/1000);
windShift = round(windShift*fs/1000);

iloc=1:windShift:length(sig)-windSize;

if (nargin==6)
    n= length(sig);
    preemphasis=1;
    sig=sig(1:n-1)-sig(2:n)*preemphasis;
    sig(n)=sig(n-1);
end


iloc=1:windShift:length(sig)-windSize;
sigbuf=[zeros(lporder,1);sig(:)];


%Append zeros to the signal to prevent the last frame from being skipped 
if iloc(end)<length(sig)
    iloc(end+1)=iloc(end)+windShift;
    length(sig)-windSize+windShift;
    sigbuf=[sigbuf(:); zeros(length(sig)-windSize+windShift,1)];
end


lpres=zeros(length(sig),1);
lpcoeff=zeros(lporder+1,length(iloc));



for i = 1:length(iloc)
    %iloc(i)
    %length(sigbuf)
    seg=sigbuf(iloc(i):iloc(i)+lporder+windSize-1);
    frmlpc=lpc(seg(lporder+1:end),lporder);
    lpcoeff(:,i)=frmlpc;
    sigest=filter([0 -frmlpc(2:end)],1,seg);
    frmlpres=filter(frmlpc,1,seg);
    ienergy(i)=sum(seg(1:end-lporder).^2)/(length(seg)-lporder);
    lpres(iloc(i):iloc(i)+windSize-1)=frmlpres(lporder+1:end);
end
lpres=lpres(1:length(sig));

    
