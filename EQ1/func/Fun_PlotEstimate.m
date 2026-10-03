function Fun_PlotEstimate(...
    G,observedMoment,m2,...
    stationLon,stationLat,...
    eventLon,eventLat,varargin)
%PLOTRESULTSFIT Plot observed and predicted second-moment results.
% INPUTS
%   G              Design matrix, N x 6.
%   observedMoment Observed temporal second moments, N x 1.
%   m2             Inverted model vector, 6 x 1:
%                  [mtt,mxt,myt,mxx,mxy,myy].
%   stationLon     Station longitudes, N x 1.
%   stationLat     Station latitudes, N x 1.
%   eventLon       Event longitude.
%   eventLat       Event latitude.
% -------------------------------------------------------------------------
% Parse optional inputs
% -------------------------------------------------------------------------

p=inputParser;

addParameter(p,'FigureNumber',[],...
    @(x) isempty(x) || (isscalar(x) && isnumeric(x)));

addParameter(p,'ColorLimits',[],...
    @(x) isempty(x) || (isnumeric(x) && numel(x)==2));

addParameter(p,'MapLimits',[],...
    @(x) isempty(x) || (isnumeric(x) && numel(x)==4));

addParameter(p,'MarkerSize',90,...
    @(x) isnumeric(x) && isscalar(x) && x>0);

addParameter(p,'ShowLabels',false,...
    @(x) islogical(x) || (isnumeric(x) && isscalar(x)));

addParameter(p,'StationNames',{},...
    @(x) iscell(x) || isstring(x));

parse(p,varargin{:});

figureNumber=p.Results.FigureNumber;
userColorLimits=p.Results.ColorLimits;
userMapLimits=p.Results.MapLimits;
markerSize=p.Results.MarkerSize;
showLabels=logical(p.Results.ShowLabels);
stationNames=p.Results.StationNames;

% -------------------------------------------------------------------------
% Format and check inputs
% -------------------------------------------------------------------------

G=double(G);
observedMoment=double(observedMoment(:));
m2=double(m2(:));
stationLon=double(stationLon(:));
stationLat=double(stationLat(:));

if size(G,2)~=6
    error('plotresultsfit:InvalidG',...
        'G must have six columns.');
end
stationCount=size(G,1);

% -------------------------------------------------------------------------
% Calculate predicted values and residuals
% -------------------------------------------------------------------------

predictedMoment=G*m2;
residualMoment=predictedMoment-observedMoment;

% Second moments should theoretically be nonnegative. Avoid complex values
% in plots if small negative values are caused by numerical roundoff.
observedMomentForPlot=max(observedMoment,0);
predictedMomentForPlot=max(predictedMoment,0);

observedTauc=2*sqrt(observedMomentForPlot);
predictedTauc=2*sqrt(predictedMomentForPlot);
residualTauc=predictedTauc-observedTauc;

% -------------------------------------------------------------------------
% Misfit statistics
% -------------------------------------------------------------------------

momentEnergy=sum(observedMoment.^2);

if momentEnergy>eps
    relativeMomentMisfit=...
        sum(residualMoment.^2)/momentEnergy;
else
    relativeMomentMisfit=NaN;
end

taucEnergy=sum(observedTauc.^2);

if taucEnergy>eps
    relativeTaucMisfit=...
        sum(residualTauc.^2)/taucEnergy;
else
    relativeTaucMisfit=NaN;
end

rmsMomentResidual=sqrt(mean(residualMoment.^2));
rmsTaucResidual=sqrt(mean(residualTauc.^2));

% -------------------------------------------------------------------------
% Color limits
% -------------------------------------------------------------------------

if isempty(userColorLimits)
    tauAll=[observedTauc;predictedTauc];
    tauAll=tauAll(isfinite(tauAll));

    if isempty(tauAll)
        colorLimits=[0,1];
    else
        colorLimits=[min(tauAll),max(tauAll)];

        if colorLimits(1)==colorLimits(2)
            colorLimits=colorLimits+[-0.05,0.05];
        end
    end
else
    colorLimits=sort(userColorLimits(:).');
end

% Residual color scale is symmetric around zero.
maxResidual=max(abs(residualTauc));

if isempty(maxResidual) || ~isfinite(maxResidual) || maxResidual==0
    residualLimits=[-1,1];
else
    residualLimits=[-maxResidual,maxResidual];
end

% -------------------------------------------------------------------------
% Map limits
% -------------------------------------------------------------------------

if isempty(userMapLimits)
    lonRange=max(stationLon)-min(stationLon);
    latRange=max(stationLat)-min(stationLat);

    if lonRange==0
        lonRange=0.1;
    end

    if latRange==0
        latRange=0.1;
    end

    paddingLon=0.08*lonRange;
    paddingLat=0.08*latRange;

    mapLimits=[...
        min(stationLon)-paddingLon,...
        max(stationLon)+paddingLon,...
        min(stationLat)-paddingLat,...
        max(stationLat)+paddingLat];

    % Ensure the event is also visible.
    mapLimits(1)=min(mapLimits(1),eventLon);
    mapLimits(2)=max(mapLimits(2),eventLon);
    mapLimits(3)=min(mapLimits(3),eventLat);
    mapLimits(4)=max(mapLimits(4),eventLat);
else
    mapLimits=userMapLimits(:).';
end

% -------------------------------------------------------------------------
% Create figure
% -------------------------------------------------------------------------

if isempty(figureNumber)
    figureHandle=figure('Color','w');
else
    figureHandle=figure(figureNumber);
    clf(figureHandle);
    set(figureHandle,'Color','w');
end

set(figureHandle,...
    'Name','Second-Moment Inversion Fit',...
    'NumberTitle','off');

% Use a 2-by-2 layout.
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% -------------------------------------------------------------------------
% Panel 1: Observed tau_c
% -------------------------------------------------------------------------
cols=slanCM('seismic',30);
nexttile;
scatter(stationLon,stationLat,...
    markerSize,observedTauc,'filled',...
    'MarkerEdgeColor','k',...
    'LineWidth',0.5);
hold on;
plot(eventLon,eventLat,'p',...
    'MarkerSize',15,...
    'MarkerFaceColor','r',...
    'MarkerEdgeColor','k');
axis(mapLimits);
colormap(gca,cols);
clim(colorLimits);
colorbar;
xlabel('Longitude');
ylabel('Latitude');
title('(a) Observed \tau_c');
Fun_defaultAxes;
% -------------------------------------------------------------------------
% Panel 2: Predicted tau_c
% -------------------------------------------------------------------------
nexttile;
scatter(stationLon,stationLat,...
    markerSize,predictedTauc,'filled',...
    'MarkerEdgeColor','k',...
    'LineWidth',0.5);
hold on;
plot(eventLon,eventLat,'p',...
    'MarkerSize',15,...
    'MarkerFaceColor','r',...
    'MarkerEdgeColor','k');
axis(mapLimits);
colormap(gca,cols);
clim(colorLimits);
colorbar;
xlabel('Longitude');
ylabel('Latitude');
title('(b) Predicted \tau_c');
Fun_defaultAxes;
% -------------------------------------------------------------------------
% Panel 3: Residual tau_c
% -------------------------------------------------------------------------
% 
% nexttile;
% 
% scatter(stationLon,stationLat,...
%     markerSize,residualTauc,'filled',...
%     'MarkerEdgeColor','k',...
%     'LineWidth',0.5);
% 
% hold on;
% plot(eventLon,eventLat,'p',...
%     'MarkerSize',15,...
%     'MarkerFaceColor','r',...
%     'MarkerEdgeColor','k');
% 
% axis equal;
% axis(mapLimits);
% grid on;
% box on;
% 
% colormap(gca,parula);
% clim(residualLimits);
% colorbar;
% 
% xlabel('Longitude');
% ylabel('Latitude');
% title('\tau_c residual: predicted - observed');
% Fun_defaultAxes;
% -------------------------------------------------------------------------
% Panel 4: Observed-predicted scatter plot
% -------------------------------------------------------------------------
nexttile;
scatter(observedTauc,predictedTauc,...
    markerSize,residualTauc,'filled',...
    'MarkerEdgeColor','k',...
    'LineWidth',0.5);
hold on;
finiteValues=isfinite(observedTauc) & isfinite(predictedTauc);
if any(finiteValues)
    minValue=min([observedTauc(finiteValues);predictedTauc(finiteValues)]);
    maxValue=max([observedTauc(finiteValues);predictedTauc(finiteValues)]);
    if minValue==maxValue
        minValue=minValue-0.05;
        maxValue=maxValue+0.05;
    end
    plot([minValue,maxValue],[minValue,maxValue],...
        'k--','LineWidth',1.2);
    xlim([minValue,maxValue]);
    ylim([minValue,maxValue]);
end
xlabel('Observed \tau_c');
ylabel('Predicted \tau_c');
title('(c) Comparision');
colormap(gca,cols);
clim(residualLimits);
colorbar;
Fun_defaultAxes;

set(gcf,'position',[300,100,1500,500])
% -------------------------------------------------------------------------
% Add station labels if requested
% -------------------------------------------------------------------------
if showLabels
    if isempty(stationNames)
        stationNames=arrayfun(@(k)sprintf('%d',k),...
            (1:stationCount).','UniformOutput',false);
    end

    if numel(stationNames)~=stationCount
        warning('plotresultsfit:StationNameCount',...
            'StationNames count does not match station count.');
    else
        % Add labels to the last axes.
        ax=gca;

        for k=1:stationCount
            text(ax,...
                observedTauc(k),...
                predictedTauc(k),...
                ['  ',char(stationNames{k})],...
                'FontSize',10,...
                'Interpreter','none');
        end
    end
end

% -------------------------------------------------------------------------
% Add summary information to figure
% -------------------------------------------------------------------------
% annotationText=sprintf([...
%     'Moment relative misfit = %.4g\n' ...
%     '\\tau_c relative misfit = %.4g\n' ],...
%     relativeMomentMisfit,...
%     relativeTaucMisfit);
% 
% annotation(figureHandle,'textbox',...
%     [0.36,0.005,0.30,0.065],...
%     'String',annotationText,...
%     'FitBoxToText','on',...
%     'BackgroundColor','w',...
%     'EdgeColor',[0.7,0.7,0.7],...
%     'FontSize',14,...
%     'HorizontalAlignment','center');

end
