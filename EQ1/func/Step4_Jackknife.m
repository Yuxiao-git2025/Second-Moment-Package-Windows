function jack=Step4_Jackknife(G,d,m2,late,lone,mlats,mlons,azband,varargin)
% Azimuthal delete-group jackknife uncertainty analysis.
%
% Inputs
% -------
% G       : data/model matrix, N-by-6
% d       : data vector, N-by-1 or 1-by-N
% m2      : best-fit model vector
% late    : event latitude
% lone    : event longitude
% mlats   : station latitudes
% mlons   : station longitudes
% azband  : azimuth bin width in degrees
%
% Optional name-value inputs
% ---------------------------
% 'MakeFigure'   : true/false, default false
% 'FigureNumber' : figure number, default 6
% 'MinObservations' : minimum observations required for inversion,
%                     default 6
%
% Output
% -------
% jack contains:
%
%   jack.m2                 All valid model realizations
%   jack.m2All              All realizations, including failed ones
%   jack.tauc
%   jack.Lc
%   jack.Wc
%   jack.Vx
%   jack.Vy
%   jack.mV0
%   jack.ratio
%   jack.usedCount
%   jack.deletedCount
%   jack.azimuthLower
%   jack.azimuthUpper
%   jack.azimuthCenter
%   jack.valid
%   jack.Nvalid
%   jack.meanModel
%   jack.covarianceSample
%   jack.covarianceJK
%   jack.errors
%
% The reported errors of derived quantities are calculated directly from
% the jackknife realizations, rather than only by linear covariance
% propagation.

% =========================================================================
% Parse optional parameters
% =========================================================================
p=inputParser;
addParameter(p,'MakeFigure',false,...
    @(x)islogical(x) || isnumeric(x));

addParameter(p,'FigureNumber',6,...
    @(x)isnumeric(x) && isscalar(x));

addParameter(p,'MinObservations',6,...
    @(x)isnumeric(x) && isscalar(x) && x>=1);

parse(p,varargin{:});

makeFigure=logical(p.Results.MakeFigure);
figureNumber=p.Results.FigureNumber;
minObservations=p.Results.MinObservations;

% =========================================================================
% Input validation
% =========================================================================

G=double(G);
d=double(d(:));
m2=double(m2(:));
mlats=double(mlats(:));
mlons=double(mlons(:));

if numel(m2)>=7
    m2=m2(1:6);
end

if size(G,1)~=numel(d)
    error('G and d have incompatible numbers of rows.');
end

if size(G,2)~=6
    error('G must have six columns.');
end

if numel(mlats)~=numel(d) || numel(mlons)~=numel(d)
    error('mlats, mlons, and d must have the same length.');
end

if numel(m2)~=6
    error('m2 must contain six model parameters.');
end

if azband<=0 || azband>360
    error('azband must be in the interval (0,360].');
end

if any(~isfinite(G(:))) || any(~isfinite(d))
    error('G and d must contain only finite values.');
end

% =========================================================================
% Basic settings
% =========================================================================
N=numel(d);
nband=ceil(360/azband);

% =========================================================================
% Calculate station azimuths
% =========================================================================
stationAzimuth=NaN(N,1);
stationDistance=NaN(N,1);
stationDelta=NaN(N,1);

for ii=1:N
    [stationDistance(ii),stationAzimuth(ii),stationDelta(ii)]=...
        distaz(late,lone,mlats(ii),mlons(ii));
end

% Convert possible 360-degree values to zero.
stationAzimuth(stationAzimuth>=360)=0;

% =========================================================================
% Preallocate
% =========================================================================

% Important:
% Use NaN, not zeros. Failed inversions must remain invalid.
m2All=NaN(nband,6);

taucAll=NaN(nband,1);
LcAll=NaN(nband,1);
WcAll=NaN(nband,1);

VxAll=NaN(nband,1);
VyAll=NaN(nband,1);
mV0All=NaN(nband,1);
ratioAll=NaN(nband,1);

usedCount=zeros(nband,1);
deletedCount=zeros(nband,1);

azimuthLower=zeros(nband,1);
azimuthUpper=zeros(nband,1);
azimuthCenter=zeros(nband,1);

% =========================================================================
% Delete-group jackknife loop
% =========================================================================

for ii=1:nband

    az1=(ii-1)*azband;
    az2=min(ii*azband,360);

    azimuthLower(ii)=az1;
    azimuthUpper(ii)=az2;
    azimuthCenter(ii)=0.5*(az1+az2);

    % Half-open intervals, except for the final interval.
    if ii<nband
        deleteIndex=...
            stationAzimuth>=az1 & stationAzimuth<az2;
    else
        deleteIndex=...
            stationAzimuth>=az1 & stationAzimuth<=az2;
    end

    useIndex=~deleteIndex;

    deletedCount(ii)=sum(deleteIndex);
    usedCount(ii)=sum(useIndex);

    if usedCount(ii)<minObservations

        warning(...
            'Jackknife bin %d skipped: only %d observations remain.',...
            ii,usedCount(ii));

        continue;
    end

    G2=G(useIndex,:);
    d2=d(useIndex);

    try

        [m2tmp,~,~,~,~,~]=Step3_Inverse(G2,d2);

        m2tmp=m2tmp(:);

        if numel(m2tmp)>=7
            m2tmp=m2tmp(1:6);
        end

        if numel(m2tmp)~=6
            error('Step3_Inverse did not return six model parameters.');
        end

        if any(~isfinite(m2tmp))
            error('Model contains NaN or Inf.');
        end

        derived=calDerived(m2tmp,G2,d2);

        % Save model parameters.
        m2All(ii,:)=m2tmp.';

        % Save derived parameters.
        taucAll(ii)=derived.tauc;
        LcAll(ii)=derived.Lc;
        WcAll(ii)=derived.Wc;

        VxAll(ii)=derived.Vx;
        VyAll(ii)=derived.Vy;
        mV0All(ii)=derived.mV0;
        ratioAll(ii)=derived.ratio;

    catch ME

        warning(...
            'Jackknife bin %d failed: %s',...
            ii,ME.message);

    end

end

% =========================================================================
% Select valid realizations
% =========================================================================

valid=all(isfinite(m2All),2);

if sum(valid)<2
    error('Fewer than two valid jackknife realizations remain.');
end

jack.m2All=m2All;
jack.valid=valid;
jack.Nvalid=sum(valid);

jack.m2=m2All(valid,:);

jack.tauc=taucAll(valid);
jack.Lc=LcAll(valid);
jack.Wc=WcAll(valid);

jack.Vx=VxAll(valid);
jack.Vy=VyAll(valid);
jack.mV0=mV0All(valid);
jack.ratio=ratioAll(valid);

jack.usedCount=usedCount;
jack.deletedCount=deletedCount;

jack.azimuthLower=azimuthLower;
jack.azimuthUpper=azimuthUpper;
jack.azimuthCenter=azimuthCenter;

jack.azimuthValid=azimuthCenter(valid);

% =========================================================================
% Model-parameter covariance
% =========================================================================

K=jack.Nvalid;

jack.meanModel=mean(jack.m2,1);

deviations=jack.m2-jack.meanModel;

% Sample covariance of jackknife realizations.
jack.covarianceSample=...
    deviations.'*deviations/(K-1);

% Classical delete-group jackknife covariance.
jack.covarianceJK=...
    (K-1)/K*(deviations.'*deviations);

% Symmetrize covariance matrices.
jack.covarianceSample=...
    0.5*(jack.covarianceSample+jack.covarianceSample.');

jack.covarianceJK=...
    0.5*(jack.covarianceJK+jack.covarianceJK.');

% =========================================================================
% Direct Jackknife uncertainty of derived quantities
% =========================================================================

jack.errors=struct();

jack.errors.sigmaTauc=jackknifeSigma(jack.tauc);
jack.errors.sigmaLc=jackknifeSigma(jack.Lc);
jack.errors.sigmaWc=jackknifeSigma(jack.Wc);

jack.errors.sigmaVx=jackknifeSigma(jack.Vx);
jack.errors.sigmaVy=jackknifeSigma(jack.Vy);
jack.errors.sigmaMv0=jackknifeSigma(jack.mV0);
jack.errors.sigmaRatio=jackknifeSigma(jack.ratio);

% Also calculate linear covariance-propagation results for comparison.
    jack.errorsLinear=getDerivedErrors(m2,jack.covarianceJK);
% =========================================================================
% Display results
% =========================================================================
file=fopen('Info.txt','a');
fprintf(file,'\n');
fprintf(file,'============================================================\n');
fprintf(file,'AZIMUTHAL JACKKNIFE RESULTS\n');
fprintf(file,'============================================================\n');
fprintf(file,'Azimuth bin             = %.1f degree\n',azband);
fprintf(file,'Total bins              = %d\n',nband);
fprintf(file,'Valid jackknife samples = %d/%d\n',jack.Nvalid,nband);
validUsedCount=usedCount(valid);
validDeletedCount=deletedCount(valid);
fprintf(file,'Minimum used arrivals   = %d\n',min(validUsedCount));
fprintf(file,'Maximum deleted arrivals= %d\n',max(validDeletedCount));

fprintf(file,'\nJackknife uncertainty:\n');
fprintf(file,'tauc = %.3f ( %.3f )\n',calDerived(m2,G,d).tauc,jack.errors.sigmaTauc);

fprintf(file,'Lc   = %.3f ( %.3f )\n',calDerived(m2,G,d).Lc,jack.errors.sigmaLc);

fprintf(file,'Wc   = %.3f ( %.3f )\n',calDerived(m2,G,d).Wc,jack.errors.sigmaWc);

fprintf(file,'Vx   = %.3f ( %.3f )\n',calDerived(m2,G,d).Vx,jack.errors.sigmaVx);

fprintf(file,'Vy   = %.3f ( %.3f )\n',calDerived(m2,G,d).Vy,jack.errors.sigmaVy);

fprintf(file,'|V0| = %.3f ( %.3f )\n',calDerived(m2,G,d).mV0,jack.errors.sigmaMv0);

fprintf(file,'L0/Lc= %.3f ( %.3f )\n',calDerived(m2,G,d).ratio,jack.errors.sigmaRatio);
fclose(file);

% =========================================================================
% Optional plots
% =========================================================================

if makeFigure

    figure(figureNumber);
    clf;
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    nexttile;
    histogram(jack.tauc,15,...
        'Normalization','pdf',...
        'FaceColor',[0.35 0.47 0.93]);
    xlabel('\tau_c');
    ylabel('PDF');
    grid on;

    nexttile;
    histogram(jack.Lc,15,...
        'Normalization','pdf',...
        'FaceColor',[0.35 0.47 0.93]);
    xlabel('L_c');
    ylabel('PDF');
    grid on;

    nexttile;
    histogram(jack.Wc,15,...
        'Normalization','pdf',...
        'FaceColor',[0.35 0.47 0.93]);
    xlabel('W_c');
    ylabel('PDF');
    grid on;

    nexttile;
    histogram(jack.Vx,15,...
        'Normalization','pdf',...
        'FaceColor',[0.35 0.47 0.93]);
    xlabel('V_x');
    ylabel('PDF');
    grid on;

    nexttile;
    histogram(jack.Vy,15,...
        'Normalization','pdf',...
        'FaceColor',[0.35 0.47 0.93]);
    xlabel('V_y');
    ylabel('PDF');
    grid on;

    nexttile;
    histogram(jack.ratio,15,...
        'Normalization','pdf',...
        'FaceColor',[0.35 0.47 0.93]);
    xlabel('L_0/L_c');
    ylabel('PDF');
    grid on;

end

end


function sigma=jackknifeSigma(values)
% Classical delete-group jackknife standard error.
values=values(:);
values=values(isfinite(values));

K=numel(values);
if K<2
    sigma=NaN;
    return;
end
meanValue=mean(values);
sigma=sqrt((K-1)/K*sum((values-meanValue).^2));
end

