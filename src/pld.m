% Estimate an apparent source-time function using projected Landweber deconvolution.
%
% The model is:
%
%     data = GF * RSTF
%
% where the convolution is performed in the FFT domain. The ASTF is
% constrained to be nonnegative and is allowed to exist only between
% samples T1 and T.
%
% The function first computes a misfit curve for a sequence of candidate
% ASTF durations. The user then selects a final duration T from the curve.
% The deconvolution is repeated using the selected T.
%
% Usage:
%     [RSTF,dhat,T,eps,t,e]=pld(data,GF,T1,niter);
% Input:
%     data    mainshock waveform, column or row vector
%     GF      empirical Green's function waveform
%     T1      first allowed ASTF sample
%     niter   number of Landweber iterations
%
% Output:
%     RSTF    estimated ASTF
%     dhat    fitted mainshock waveform, GF * RSTF
%     T       selected ASTF end sample
%     eps     normalized final waveform misfit
%     t       candidate ASTF end samples
%     e       normalized misfit for each candidate T
function [RSTF,dhat,T,fiterror,t,e]=pld(data,GF,T1,niter)
data=data(:);
GF=GF(:);
if nargin<4 || isempty(niter)
    niter=50;
end
niter=round(niter);
N=max(numel(data),numel(GF));
% Check zero-dimension
if numel(data)<N
    data(end+1:N)=0;
end
if numel(GF)<N
    GF(end+1:N)=0;
end
T1=round(T1);
T1=max(1,min(T1,N));
Energy=sum(data.^2);
GFw=fft(GF);
dataw=fft(data);
spectralNorm=max(abs(GFw).^2);

if ~isfinite(spectralNorm) || spectralNorm<=eps
    error('pld:ZeroGF','The EGF has zero or negligible spectral energy.');
end

% Landweber step length
step=1/spectralNorm;
% Candidate ASTF end samples
% The original implementation searched from T1/N to 0.4 of the record
minimumFraction=T1/N;
maximumFraction=0.4;
fractionStep=0.005;
if minimumFraction>maximumFraction
    error('pld:InvalidDurationRange','T1 is too late for the selected FFT record length.');
end
candidateFraction=minimumFraction:fractionStep:maximumFraction;
candidateT=round(candidateFraction*N);
candidateT=max(candidateT,T1);
candidateT=min(candidateT,N);
candidateT=unique(candidateT);

t=candidateT(:);
e=nan(numel(t),1);

% Compute the misfit curve
for k=1:numel(t)
    [~,dhatTrial]=Fun_SolvePLD(data,dataw,GFw,T1,t(k),niter,step);

    residual=data-dhatTrial;
    e(k)=sum(residual.^2)/Energy;
end
% Display the duration-misfit curve and select the final duration
figure('Name','PLD duration-misfit curve','NumberTitle','off');
nexttile;
plot(t,e,'k','LineWidth',1.2);
hold on;
xlabel('ASTF end sample');
ylabel('Normalized misfit');
title('Select the ASTF duration');
Fun_defaultAxes;
xlim([min(t),max(t)]);
ylim([0,max(1,max(e(isfinite(e))))]);

[xPick,~]=ginput(1);
if isempty(xPick) || ~isfinite(xPick)
    close(gcf);
    error('pld:Cancelled','No ASTF duration was selected.');
end
[~,selectedIndex]=min(abs(t-xPick));
T=t(selectedIndex);
xline(T,'Color',[.65 .65 .65],'LineStyle',':','LineWidth',1.4);
drawnow;


% Recompute the final ASTF with the selected duration
[RSTF,dhat]=Fun_SolvePLD(data,dataw,GFw,T1,T,niter,step);
residual=data-dhat;
fiterror=sum(residual.^2)/Energy;
if isgraphics(gcf)
    close(gcf);
end

end


% Solve the projected Landweber problem for a fixed ASTF end sample
function [RSTF,dhat]=Fun_SolvePLD(data,dataw,GFw,T1,T,niter,step)

N=numel(data);

T=round(T);
T=max(T1,min(T,N));

source=zeros(N,1);
sourcew=fft(source);

GFconj=conj(GFw);

for iteration=1:niter
    modelw=GFw.*sourcew;
    residualw=dataw-modelw;
    gradientw=GFconj.*residualw;
    gradient=real(ifft(gradientw));
    source=source+step*gradient;
    % Projection onto:
    %   1) the nonnegative cone;
    %   2) the allowed time interval [T1,T]
    source(1:T1-1)=0;
    source(T+1:N)=0;
    source=max(source,0);
    sourcew=fft(source);
end
RSTF=source;
dhat=real(ifft(GFw.*sourcew));

end



% function[RSTF,dhat,T,eps,t,e]= pld(data,GF,T1,niter)
% % DOES THE PLD inversion process for the RSTF from a data seismogram
% % and a GF seismogram.   The P or S-wave seismograms are expected to be
% % alligned so that the P-wave arrival time occurs in the sample
% % T1, and the GF begins at the P-wave arrival time.
% % Also, it is assumed that they have been padded with zeros to
% % a sufficient length of time to avoid wrap-arounds. (i.e. more than
% % double).  T-T1 is the estimate of the duration of the STF
% % you pick fromt the misfit tradeoff curve.
% % niter is the number of iterations (ie inverse damping)
% %
% %   The returned variables are:
% %   RSTF  the apparent source time function
% %   dhat  the fit to the data seismogram (RSTF*GF)
% %   T     the duration pick
% %   eps   the misfit at T
% %   t,e   the misfit vs duration tradeoff curve
% %
% %  This follows Bertero et al. 1997, Lanza et al. 1999, and McGuire 2004,
% %  see manual for references.
% 
% N=length(data);
% e0=(T1)/N;
% epsilon=e0:.005:0.4;
% ne=length(epsilon);
% e=zeros(ne,1);
% t=zeros(ne,1);
% foldt=zeros(N,1);
% fnewt=zeros(N,1);
% fneww=fft(fnewt);
% GT=zeros(N,1);
% dataw=fft(data);
% GFw=fft(GF);
% Gstarw=conj(GFw);
% tau=max(abs(GFw));
% tau=tau^2;
% tau=1/tau;
% 
% for i=1:N
%     GT(i)=GF(N-i+1);
% end
% GTw=fft(GT);
% nit=niter;
% % LOOP OVER EPSILONS
% 
% for i=1:ne
%     eps=epsilon(i);
%     T=round(eps*N);
%     t(i)=T;
%     %set up inverse problem
% 
%     % iterate over n
%     for j=1:nit
%         foldt=fnewt;
%         foldw=fneww;
%         clear dum
%         f=fft(foldt);
%         res=dataw-GFw.*f;
%         dum=Gstarw.*res;
%         gneww=foldw+tau.*dum;
%         gnewt=real(ifft(gneww));
%         fnewt=posproj(gnewt,T1,T);
%         fneww=fft(fnewt);
%     end
%     % end iterations for this T.
%     dhat=real(ifft(fneww.*GFw));
% 
%     sum1=0; sum2=0;
%     for j=1:N
%         sum1=sum1+(data(j)-dhat(j))^2;
%         sum2=sum2+(data(j))^2;
%     end
%     e(i)=sum1/sum2;
% end
% 
% % plot eps vs T
% 
% figure;
% plot(t,e,'Color','k','LineWidth',0.8);
% hold on;
% ylabel('Epsilon')
% xlabel('T');
% axis([0 round(epsilon(end)*N) 0 1])
% [x,~] = ginput(1);
% T=round(x);
% xline(T,'LineWidth',1.6,'Color',[0.65 0.65 0.65],'LineStyle',':');
% % iterate over n
% foldt=zeros(N,1);
% fnewt=zeros(N,1);
% fneww=fft(fnewt);
% for j=1:nit
%     foldt=fnewt;
%     foldw=fneww;
%     clear dum
%     f=fft(foldt);
%     res=dataw-GFw.*f;
%     dum=Gstarw.*res;
%     gneww=foldw+tau.*dum;
%     gnewt=real(ifft(gneww));
%     fnewt=posproj(gnewt,T1,T);
%     fneww=fft(fnewt);
% end
% 
% dhat=real(ifft(fneww.*GFw));
% RSTF=fnewt;
% sum1=0; sum2=0;
% for j=1:N
%     sum1=sum1+(data(j)-dhat(j))^2;
%     sum2=sum2+(data(j))^2;
% end
% eps=sum1/sum2;
% close
