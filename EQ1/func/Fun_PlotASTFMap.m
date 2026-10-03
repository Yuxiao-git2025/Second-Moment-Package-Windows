function hFig=Fun_PlotASTFMap(IJ,slon,slat,t2,STF,dtsv,t0,lone,late)
% Plot station characteristic durations on a longitude-latitude map
% and overlay a small ASTF waveform near each station.
%
% Inputs:
%   IJ      station indices to plot
%   slon    station longitudes
%   slat    station latitudes
%   t2      second central moment, in s^2
%   STF     ASTF matrix, one station per row
%   dtsv    sampling interval for each station, in s
%   t0      ASTF centroid position in MATLAB sample-index coordinates
%   lone    longitude vector for an additional line
%   late    latitude vector for an additional line
%
% Important:
%   Because t0 is based on MATLAB indices starting from 1,
%   its physical time is:
%
%       t0Time=(t0-1)*dt
%
%   The input t2 must already be converted to s^2.
% =========================
% Check input dimensions
% =========================
if nargin<9
    error('Fun_PlotASTFMap requires 9 input arguments.');
end
IJ=IJ(:);
slon=slon(:);
slat=slat(:);
t2=t2(:);
dtsv=dtsv(:);
t0=t0(:);
if any(~isfinite(IJ)) || any(IJ~=round(IJ))
    error('IJ must contain finite integer indices.');
end
IJ=round(IJ);
nStation=numel(slon);

if numel(slat)~=nStation
    error('slon and slat must have the same length.');
end

if numel(t2)~=nStation || ...
        numel(dtsv)~=nStation || ...
        numel(t0)~=nStation
    error('Station-dependent inputs have incompatible lengths.');
end
% Remove invalid station indices from plotting.
validStation=isfinite(slon(IJ)) & ...
    isfinite(slat(IJ)) & ...
    isfinite(t2(IJ)) & ...
    t2(IJ)>=0 & ...
    isfinite(dtsv(IJ)) & ...
    dtsv(IJ)>0;

if ~any(validStation)
    error('No valid station has finite coordinates and t2.');
end

IJ=IJ(validStation);

% =========================
% Characteristic duration
% =========================

tauc=nan(size(t2));
validT2=isfinite(t2) & t2>=0;
tauc(validT2)=2*sqrt(t2(validT2));

validDuration=isfinite(tauc(IJ)) & tauc(IJ)>=0;

if ~any(validDuration)
    error('No valid characteristic duration is available.');
end

taucMedian=median(tauc(IJ(validDuration)),'omitnan');

if ~isfinite(taucMedian) || taucMedian<=0
    taucMedian=1;
end

% =========================
% Create the map
% =========================

hFig=figure(...
    'Name','ASTF spatial distribution',...
    'NumberTitle','off',...
    'Color','w');

mapAx=axes('Parent',hFig);
hold(mapAx,'on');

durationValues=tauc(IJ);

hScatter=scatter(mapAx,slon(IJ),slat(IJ),100,durationValues,'filled','MarkerEdgeColor',[.2 .2 .2]);

% Keep the scatter object in case the caller wants to modify it.
setappdata(hFig,'DurationScatter',hScatter);
% Colormap.
colormap(mapAx,slanCM('seismic',30));
hColorbar=colorbar(mapAx);
hColorbar.Label.String='Characteristic Duration (s)';

finiteColor=durationValues(isfinite(durationValues));
if isempty(finiteColor)
    colorMin=0;
    colorMax=1;
else
    colorMin=min(finiteColor);
    colorMax=max(finiteColor);
    if colorMin==colorMax
        colorMin=colorMin-0.5;
        colorMax=colorMax+0.5;
    else
        colorRange=colorMax-colorMin;
        colorMin=colorMin-0.1*colorRange;
        colorMax=colorMax+0.1*colorRange;
    end
end
clim(mapAx,[colorMin,colorMax]);

% Plot the additional line or boundary.
if ~isempty(lone) && ~isempty(late)
    if numel(lone)~=numel(late)
        error('lone and late must have the same number of elements.');
    end

    plot(mapAx,lone,late,...
        'k-',...
        'LineWidth',1.0,...
        'Marker','pentagram',...
        'MarkerSize',12,...
        'MarkerFaceColor','r');
end

if exist('Fun_Decorat','file')==2
    Fun_Decorat;
end

if exist('Fun_SetTickNumber','file')==2
    Fun_SetTickNumber(mapAx,4);
end

box(mapAx,'off');

set(hFig,...
    'Position',[400,85,800,680]);

drawnow;

% =========================
% Convert station coordinates
% to normalized figure units
% =========================

nPlot=numel(IJ);
Xf=nan(nPlot,1);
Yf=nan(nPlot,1);

for k=1:nPlot
    [Xf(k),Yf(k)]=ds2nfu(mapAx,...
        slon(IJ(k)),...
        slat(IJ(k)));
end

% Size of the small ASTF axes in normalized figure units.
smallWidth=0.08;
smallHeight=0.13;

Npts=size(STF,2);

% =========================
% Overlay ASTF waveforms
% =========================

for k=1:nPlot

    stationIndex=IJ(k);

    % Skip stations with invalid ASTF data.
    if stationIndex<1 || stationIndex>size(STF,1)
        continue;
    end

    dt=dtsv(stationIndex);
    astf=STF(stationIndex,:);

    if ~isfinite(dt) || dt<=0 || isempty(astf)
        continue;
    end

    % Physical time axis.
    tAxis=(0:Npts-1)*dt;

    % Get the station ASTF.
    astf=astf(:).';

    % Ensure consistency with Npts.
    if numel(astf)<Npts
        astf(end+1:Npts)=0;
    elseif numel(astf)>Npts
        astf=astf(1:Npts);
    end

    % Convert centroid from sample index to physical time.
    if isfinite(t0(stationIndex))
        t0Time=(t0(stationIndex)-1)*dt;
    else
        t0Time=NaN;
    end

    % Select a local plotting window.
    if isfinite(t0Time) && isfinite(taucMedian)
        x1=t0Time-2*taucMedian;
        x2=t0Time+3*taucMedian;

        if x2<=x1
            x1=0;
            x2=max(tAxis);
        end
    else
        x1=0;
        x2=max(tAxis);
    end

    % If x2 is not useful, use the complete ASTF window.
    if ~isfinite(x1) || ~isfinite(x2) || x2<=x1
        x1=0;
        x2=max(tAxis);
    end

    % Keep the small axes inside the figure.
    xPos=min(max(Xf(k),0),1-smallWidth);
    yPos=min(max(Yf(k),0),1-smallHeight);

    astfAx=axes(...
        'Parent',hFig,...
        'Units','normalized',...
        'Position',[xPos,yPos,smallWidth,smallHeight],...
        'Color','none',...
        'Box','off',...
        'XColor','none',...
        'YColor','none',...
        'HitTest','off',...
        'HandleVisibility','off');

    hold(astfAx,'on');

    % Replace invalid ASTF samples for plotting only.
    astfPlot=astf;
    astfPlot(~isfinite(astfPlot))=0;

    plot(astfAx,tAxis,astfPlot,...
        'k-',...
        'LineWidth',0.6);

    xlim(astfAx,[x1,x2]);

    % Set a stable y-axis range for zero or nearly zero ASTFs.
    finiteASTF=astfPlot(isfinite(astfPlot));

    if isempty(finiteASTF) || max(abs(finiteASTF))<=eps
        ylim(astfAx,[-1,1]);
    else
        yMax=max(abs(finiteASTF));

        if yMax<=0 || ~isfinite(yMax)
            yMax=1;
        end

        ylim(astfAx,[-0.05*yMax,1.05*yMax]);
    end

    axis(astfAx,'off');
end

% Make sure the map axes stay at the bottom layer.
uistack(mapAx,'bottom');

end
