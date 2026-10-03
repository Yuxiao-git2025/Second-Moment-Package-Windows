% Calculate the temporal second central moment of an apparent source-time function.
%
% The input STF is treated as a nonnegative moment-rate-like function.
% The output t2 is calculated in sample^2. The caller should multiply t2
% by dt^2 to convert it to s^2.
%
% The output t0 is the temporal centroid in sample index coordinates.
% Since the centroid is a weighted average, t0 is generally not an integer.
% Output:
%     t2      temporal second central moment in sample^2
%     t1      legacy output, equal to t0
%     t0      temporal centroid in sample index coordinates
%
function [t2,t1,t0]=findt2(STF,plott2)
STF=STF(:);
if nargin<2 || isempty(plott2)
    plott2=0;
end
N=numel(STF);
if plott2==1
    figure('Name','Select ASTF integration interval','NumberTitle','off');
    plot(1:N,STF,'k','LineWidth',1.2);
    hold on;
    xlabel('Sample');
    ylabel('ASTF amplitude');
    title('Select the ASTF integration interval');
    Fun_defaultAxes;
    xlim([1,N]);

    [xLeft,~]=ginput(1);
    if isempty(xLeft) || ~isfinite(xLeft)
        close(gcf);
        error('findt2:Cancelled','The left integration limit was not selected.');
    end

    LI=round(xLeft);
    LI=max(1,min(LI,N));
    xline(LI,'Color',[.65 .65 .65],'LineStyle',':','LineWidth',1.4);

    [xRight,~]=ginput(1);
    if isempty(xRight) || ~isfinite(xRight)
        close(gcf);
        error('findt2:Cancelled','The right integration limit was not selected.');
    end

    RI=round(xRight);
    RI=max(1,min(RI,N));
    xline(RI,'Color',[.65 .65 .65],'LineStyle',':','LineWidth',1.4);

    IntRange=sort([LI,RI]);
    I1=IntRange(1);
    I3=IntRange(2);
else
    I1=1;
    I3=N;
end

sampleIndex=(I1:I3).';
weight=STF(I1:I3);
momentSum=sum(weight);

if ~isfinite(momentSum) || momentSum<=0
    if plott2==1 && exist('gcf','var') && isgraphics(gcf)
        close(gcf);
    end
    error('findt2:InvalidMoment','The integrated ASTF moment must be positive.');
end

% Temporal centroid in sample-index coordinates.
t0=sum(sampleIndex.*weight)/momentSum;

if ~isfinite(t0)
    if plott2==1 && isgraphics(gcf)
        close(gcf);
    end
    error('findt2:InvalidCentroid','The ASTF temporal centroid is invalid.');
end

% Second central moment in sample^2.
t2=sum(weight.*(sampleIndex-t0).^2)/momentSum;

if ~isfinite(t2) || t2<0
    if plott2==1 && isgraphics(gcf)
        close(gcf);
    end
    error('findt2:InvalidSecondMoment','The ASTF second central moment is invalid.');
end

t1=t0;
if plott2==1
    t0Index=round(t0);
    t0Index=max(1,min(t0Index,N));

    plot(t0,STF(t0Index),'p',...
        'MarkerFaceColor',[1 .1 .1],...
        'MarkerEdgeColor',[1 .1 .1],...
        'MarkerSize',18);
    title(['Selected interval, tau = ',num2str(2*sqrt(t2),'%.4g'),' samples']);
    drawnow;
    pause(0.3);
    if isgraphics(gcf)
        close(gcf);
    end
end

end



% function [SecondMoment,t1sample,t0sample]=findt2(STF,plott2)
% %finds the 2nd moment of a RSTF (t2) and mean/centroid time (t1)
% % INPUTS
% % STF is the source time function vector.
% %
% % IF plott2=1, it will plot the time function and let you pick the start
% % and endpoints to integrate over.   
% % IF plott2=0 it assumes the whole time function is relaible and integrates
% % from start to end of the array.
% if(plott2)
%  figure;hold on;
%  plot(STF,'Color','k','LineWidth',0.8);
%  axis([0 round(.5*length(STF)) 0 max(STF)]);
%  disp('Pick lefthand point to start variance calculation at')
%  [x, y] = ginput(1);
%  I1=round(x);
%  xline(I1,'LineWidth',1.6,'Color',[0.65 0.65 0.65],'LineStyle',':');
%  disp('Pick righthand point to end variance calculation at')
%  [x, y] = ginput(1);
%  I3=round(x);
%  xline(I3,'LineWidth',1.6,'Color',[0.65 0.65 0.65],'LineStyle',':');
% else
%  I1=1;
%  I3=length(STF);
% end
% 
% 
% % moment and centroid
% msum=0; 
% tsum=0;
% for i=I1:I3
%   msum=msum+STF(i);
%   tsum=tsum+STF(i)*i;
% end
% % Time
% t0sample=tsum/msum;
% t1sample=tsum/msum;
% tsum=0;
% for i=I1:I3
%   tsum=tsum+STF(i)*(i-t0sample)^2; % unit: (timeindex)^2
% end 
% % Moment
% SecondMoment=tsum/msum;
% 
% if(plott2)
%  hold on;
%  plot(t0sample,STF(round(t0sample)),'Marker','pentagram','MarkerFaceColor',[1 .1 .1],'MarkerSize',20);
%  axis([0 700 0 max(STF)])
%  pause(0.2);
%  close
% end
% 
% return
% end
