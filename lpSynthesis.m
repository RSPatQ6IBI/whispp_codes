function[sig]=lpSynthesis(lpcoeff,lpres,fs,windSize,windShift)

%
%function[sig]=lpSynthesis(lpcoeff,lpres,fs,windSize,windShift,preemphasis)
%
% Synthesizes the signal given the residual and lpcoefficients. If the
% original residual and lpc's are provided, the signal should be
% recoverable without loss except for finite precision errors which are
% significantly smaller than that obtained through 24 bit quantization.
%
% This function works best with the parameters returned by lpAnalysis
%

%
% Coded by: Anand Joseph, SVL, IIIT Hyderabad
% Revision: 1a, 15/07/2009
%

windSize = round(windSize*fs/1000);
windShift = round(windShift*fs/1000);

iloc=1:windShift:length(lpres)-windSize;

%Append zeros to the signal to prevent the last frame from being skipped 
if iloc(end)<length(lpres)
    iloc(end+1)=iloc(end)+windShift;
    lpres2=[lpres(:); zeros(length(lpres)-windSize+windShift,1)];
end

sig=zeros(iloc(end)+windSize,1);


%%%%% Determine no of coeffs in lpcoeff %%%%%%
[noCoeffs,noSegs]=size(lpcoeff);
lporder=min(noCoeffs,noSegs);

% noSegs=
% if noSegs~=length(iloc)
%     lporder=noSegs-1;
%     noSegs=noCoeffs;
%     noCoeffs=lporder+1;
% else
%     lporder=noCoeffs-1;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


prevresframe=0;
iloc=iloc;
% length(iloc)
for i=1:length(iloc)-2;
    resframe=[prevresframe; lpres2(iloc(i):iloc(i+1)-1)];
    seg_est=filter(1,lpcoeff(:,i),resframe);
    seg_est=seg_est(1+length(prevresframe):end);
    sig(iloc(i):iloc(i+1)-1)=seg_est;
    prevresframe=filter(lpcoeff(:,i+1),1,seg_est);
end
% length(lpres)
% length(sig)
sig=sig(1:length(lpres));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reason for computing prevresframe using the LPC's of the next frame %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In LP analysis,  the current sample is computed as the weighted sum %
% of the  previous  p  samples. This  means  that  for  computing the %
% current sample, we need the previous p samples.  Since the previous %
% p samples have to be passed to the synthesis filter, we are passing %
% the currently estimated signal segment  through an  analysis filter %
% whose coefficients are the LPCs computed for the next frame. In the %
% next  iteration,  these  will be  passed  back to the LPCs and  the %
% original signal of the previous frame will be reconstructed.        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
