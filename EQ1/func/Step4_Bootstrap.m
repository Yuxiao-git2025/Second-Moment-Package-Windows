function [mv0u,mv0l,bound2u,bound2l,Lcu,Lcl,taucu,taucl,boot]=Step4_Bootstrap(A,b,conf,NB,varargin)
% Bootstrap uncertainty analysis for the second-moment inversion.
%
% Inputs
% -------
% A    : Green's-function / partial-derivative matrix, N-by-6
% b    : data vector, N-by-1 or 1-by-N
% conf : confidence level, e.g. 0.95
% NB   : number of bootstrap realizations
%
% Optional name-value inputs
% ---------------------------
% 'MakeFigure'  : true or false, default true
% 'FigureNumber': figure number, default 11
% 'RandomSeed'  : random seed, default []
%
% Outputs
% -------
% The first 11 outputs preserve the original bootstrap2nds format.
%
% boot contains:
%   boot.m2
%   boot.tauc
%   boot.Lc
%   boot.Wc
%   boot.Vx
%   boot.Vy
%   boot.mV0
%   boot.Vc
%   boot.L0
%   boot.ratio
%   boot.bound2
%   boot.minvr
%   boot.valid
%   boot.Nvalid
%   boot.confidenceIntervals

% -------------------------------------------------------------------------
% Parse optional parameters
% -------------------------------------------------------------------------

p=inputParser;

addParameter(p,'MakeFigure',true,@(x)islogical(x) || isnumeric(x));
addParameter(p,'FigureNumber',11,@(x)isnumeric(x) && isscalar(x));
addParameter(p,'RandomSeed',[],@(x)isempty(x) || ...
    (isnumeric(x) && isscalar(x)));

parse(p,varargin{:});

makeFigure=logical(p.Results.MakeFigure);
figureNumber=p.Results.FigureNumber;
randomSeed=p.Results.RandomSeed;

if ~isempty(randomSeed)
    rng(randomSeed);
end

% -------------------------------------------------------------------------
% Input validation
% -------------------------------------------------------------------------

A=double(A);
b=double(b(:));

if size(A,1)~=numel(b)
    error('A and b have incompatible numbers of rows.');
end

if size(A,2)~=6
    error('A must have six columns.');
end

if NB<1 || NB~=round(NB)
    error('NB must be a positive integer.');
end

if conf<=0 || conf>=1
    error('conf must be between 0 and 1.');
end
N=numel(b);
% -------------------------------------------------------------------------
% Preallocate bootstrap output
% -------------------------------------------------------------------------

boot.m2=NaN(NB,6);

boot.tauc=NaN(NB,1);
boot.Lc=NaN(NB,1);
boot.Wc=NaN(NB,1);

boot.Vx=NaN(NB,1);
boot.Vy=NaN(NB,1);
boot.mV0=NaN(NB,1);

boot.Vc=NaN(NB,1);
boot.L0=NaN(NB,1);
boot.ratio=NaN(NB,1);

boot.bound2=NaN(NB,1);
boot.minvr=NaN(NB,1);

boot.Misfit=NaN(NB,1);

% Generate all resampling indices at once.
sampleIndex=randi(N,NB,N);

% -------------------------------------------------------------------------
% Bootstrap loop
% -------------------------------------------------------------------------
for ii=1:NB
    index=sampleIndex(ii,:);
    Gi=A(index,:);
    di=b(index);
    try

        % If Step3_Inverse accepts column vectors, use:
        % m2tmp=Step3_Inverse(Gi,di);
        %
        % The transpose below preserves compatibility with your original
        % code, where the data vector was passed as a row vector.
        m2tmp=Step3_Inverse(Gi,di.');

        m2tmp=m2tmp(:);

        % Remove possible dummy variable.
        if numel(m2tmp)>=7
            m2tmp=m2tmp(1:6);
        end

        if numel(m2tmp)~=6
            error('Step3_Inverse did not return six model parameters.');
        end

        % Calculate all derived quantities.
        derived=calDerived(m2tmp,Gi,di);

        % Save model parameters.
        boot.m2(ii,:)=m2tmp.';

        % Save derived quantities.
        boot.tauc(ii)=derived.tauc;
        boot.Lc(ii)=derived.Lc;
        boot.Wc(ii)=derived.Wc;

        boot.Vx(ii)=derived.Vx;
        boot.Vy(ii)=derived.Vy;
        boot.mV0(ii)=derived.mV0;

        boot.Vc(ii)=derived.Vc;
        boot.L0(ii)=derived.L0;
        boot.ratio(ii)=derived.ratio;

        % This is the quantity used in the original bootstrap2nds.
        boot.bound2(ii)=0.5*derived.Vc;
        % Original definition:
        % minvr=max(bound2,mv0)
        boot.minvr(ii)=max(boot.bound2(ii),boot.mV0(ii));
        boot.Misfit(ii)=derived.Misfit;
    catch ME
        % Keep this realization as NaN and continue.
        fprintf('Bootstrap realization %d failed: %s\n',ii,ME.message);

    end
end

% -------------------------------------------------------------------------
% Select valid realizations
% -------------------------------------------------------------------------

boot.valid=all(isfinite(boot.m2),2);
boot.Nvalid=sum(boot.valid);

if boot.Nvalid<2
    error('Fewer than two valid bootstrap realizations remain.');
end

% Remove invalid rows from all stored quantities.
fieldNames={...
    'm2','tauc','Lc','Wc','Vx','Vy','mV0',...
    'Vc','L0','ratio','bound2','minvr','Misfit'};

for ii=1:numel(fieldNames)
    fieldName=fieldNames{ii};
    boot.(fieldName)=boot.(fieldName)(boot.valid,:);
end

% -------------------------------------------------------------------------
% Percentile confidence intervals
% -------------------------------------------------------------------------

alpha=(1-conf)/2;

probabilities=[alpha,1-alpha];

boot.confidenceLevel=conf;
boot.confidenceProbabilities=probabilities;

boot.taucCI=quantile(boot.tauc,probabilities);
boot.LcCI=quantile(boot.Lc,probabilities);
boot.WcCI=quantile(boot.Wc,probabilities);

boot.VxCI=quantile(boot.Vx,probabilities);
boot.VyCI=quantile(boot.Vy,probabilities);
boot.mV0CI=quantile(boot.mV0,probabilities);

boot.VcCI=quantile(boot.Vc,probabilities);
boot.L0CI=quantile(boot.L0,probabilities);
boot.ratioCI=quantile(boot.ratio,probabilities);

boot.bound2CI=quantile(boot.bound2,probabilities);
boot.minvrCI=quantile(boot.minvr,probabilities);

boot.MisfitCI=quantile(boot.Misfit,probabilities);

% -------------------------------------------------------------------------
% Preserve the original output convention
% -------------------------------------------------------------------------

mv0l=boot.mV0CI(1);
mv0u=boot.mV0CI(2);

Lcl=boot.LcCI(1);
Lcu=boot.LcCI(2);

taucl=boot.taucCI(1);
taucu=boot.taucCI(2);

bound2l=boot.bound2CI(1);
bound2u=boot.bound2CI(2);

minvrl=boot.minvrCI(1);
minvru=boot.minvrCI(2);

% -------------------------------------------------------------------------
% Print results
% -------------------------------------------------------------------------
clc;
file=fopen('Info.txt','w'); % Renewed every time
fprintf(file,'\n');
fprintf(file,'============================================================\n');
fprintf(file,'BOOTSTRAP CONFIDENCE INTERVALS\n');
fprintf(file,'============================================================\n');
fprintf(file,'Confidence level: %.2f %%\n',100*conf);
fprintf(file,'Valid samples   : %d/%d\n\n',boot.Nvalid,NB);

printCI('tauc',boot.taucCI,file);
printCI('Lc',boot.LcCI,file);
printCI('Wc',boot.WcCI,file);
printCI('Vx',boot.VxCI,file);
printCI('Vy',boot.Vy,file);
printCI('mV0',boot.mV0CI,file);
printCI('Vc',boot.VcCI,file);
printCI('L0',boot.L0CI,file);
printCI('L0/Lc',boot.ratioCI,file);
% printCI('bound2',boot.bound2CI,file);
% printCI('minvr',boot.minvrCI,file);
fclose(file);
% -------------------------------------------------------------------------
% Plot bootstrap distributions
% -------------------------------------------------------------------------

if makeFigure

    figure(figureNumber);
    clf;
    nums=12;
    arr=tiledlayout(2,3,"TileSpacing",'compact','Padding','compact');
    ylabel(arr,'PDF','FontSize',18,'FontName','Times New Roman');
    nexttile;
    histogram(boot.Lc,nums,'FaceColor',[0.35 0.47 0.93],'Normalization','pdf','LineWidth',0.5);
    xlabel('$L_c$','Interpreter','latex');
    grid on;Fun_defaultAxes;

    nexttile;
    histogram(boot.Wc,nums,'FaceColor',[0.35 0.47 0.93],'Normalization','pdf','LineWidth',0.5);
    xlabel('$W_c$','Interpreter','latex');
    grid on;Fun_defaultAxes;

    nexttile;
    histogram(boot.tauc,nums,'FaceColor',[0.35 0.47 0.93],'Normalization','pdf','LineWidth',0.5);
    xlabel('$\tau_c$','Interpreter','latex');
    grid on;Fun_defaultAxes;

    nexttile;
    histogram(boot.mV0,nums,'FaceColor',[0.35 0.47 0.93],'Normalization','pdf','LineWidth',0.5);
    xlabel('$|V_0|$','Interpreter','latex');
    grid on;Fun_defaultAxes;

    nexttile;
    histogram(boot.ratio,nums,'FaceColor',[0.35 0.47 0.93],'Normalization','pdf','LineWidth',0.5);
    xlabel('$L_0/L_c$','Interpreter','latex');
    grid on;Fun_defaultAxes;

    nexttile;
    histogram(boot.Misfit,nums,'FaceColor',[0.35 0.47 0.93],'Normalization','pdf','LineWidth',0.5);
    xlabel('Misfit');
    grid on;Fun_defaultAxes;

    set(gcf,'position',[300,200,1200,600]);
end

end

function printCI(name,limits,file)
fprintf(file,'%-6s : [%.2f, %.2f]\n',name,limits(1),limits(2));

end
