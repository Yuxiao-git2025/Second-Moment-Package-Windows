% >> Control parameters
fprintf('[Now] Doing Time: %s \n',datetime('now'));
fclose all;
isData=1;       % 1 to just load the seismograms from an example .mat file (*default*)
                % 0 to read in miniseed files
isMeasure=0;    % 1 to go through the EGF deconvolutions station by station;
                % 0 just loads example measurements .mat file
isInt=0;        % 0 calculate mu^(0,2) (s) by integrating over whole ASTF (*default*)
                % 1 if you want to have add a second round of interactive
                % picking of the integration limits.
isplot=1;       % 1 to make plots of the resulting STFs and data fits
isInverse=1;    % 1 to run the convex optimization inversion (*default*)
niter=80;       % Iteration numbers (20-100 is Okay)


% >> DATA INPUT
% Do your data input here;
% stasm   list of stations to try measurements on
% compm   associated list of individual components
% velMS   mainshock velocity seismograms (same order as stasm compm)
% velEGF  EGF velocity seismograms (same order as stasm compm)
% slat,slon   station latitude+longitudes
% late,lone,depe event location
if isData
    load('Input\Data.mat');
else
    Step1_Loading;
end
fprintf('[1] DONE WITH DATA-PREPARE\n');

% Parameters for the PLD measurements of ASTF duration
% t2: the temporal second moment (i.e., mu(0,2))
% Done: whether the calculation for this channel has been completed
% STF: obtained by deconvolving the mainshock waveform and the EGF
% GFsv: EGF waveform after truncation, alignment, and zero-padding
% dhatsv: fitting waveform by PLD
% datasv: MS waveform after truncation and other processing
% Tsv: duration time when PLD and T1sv is reference arrival time
% epsv: final misfit when PLD
% epldsv: misfit curves when PLD (changes with ASTF duration)
% tpldsv: duration curves when PLD (see epldsv)
% t0=t1: centroid time in ASTF (index)
% Phasesv: phase info selected (P/S)
if isMeasure
    [t2,DONE,STF,GFsv,dhatsv,datasv,Tsv,T1sv,epsv,epldsv,tpldsv,t0,t1,PhaseSv]=Step2_Measure( ...
        velEGF,velMS,niter,npMS,npEGF,dtsv,stasm,compm,isInt);
else
    load('Input\MEASUREMENTS_Exmaple.mat');
end
IJ=find(DONE==1); % effective index
% Plot the station summary
if isplot
    Fun_PlotASTFMap(IJ,slon,slat,t2,STF,dtsv,t0,lone,late);
end
fprintf('[2] DONE WITH MEASUREMENTS\n');


%%
% Aming function: d=Gm
% Read station locations, second moments, and P/S phase data
%   ↓
% Read velocity model Vp, Vp/Vs
%   ↓
% Compute local coordinates of the source and stations
%   ↓
% Calculate distance and azimuth from source to each station
%   ↓
% Call the TOPP ray tracing fortran program
%   ↓
% Obtain ray takeoff angle and travel time
%   ↓
% Construct ray direction and slowness vector
%   ↓
% Rotate slowness into fault coordinate system
%   ↓
% Construct G matrix
%   ↓
% Solve the semi-definite constrained least squares problem
%   ↓
% Obtain m2
%   ↓
% Calculate tau_c, L_c, W_c, v0, L0/Lc
clc;
if isInverse
    % =========================
    % Assemble station data
    % =========================
    IJ=IJ(:);
    Nj=numel(IJ);
    mlats=slat(IJ);
    mlons=slon(IJ);
    melevs=zeros(Nj,1);
    d=t2(IJ);
    phas=upper(char(PhaseSv(IJ)));
    % Force to horizontal quantity
    mlats=mlats(:)';
    mlons=mlons(:)';
    melevs=melevs(:)';
    d=d(:)'; % Observations
    phas=phas(:)'; % Selected phases
    % =========================
    % Load velocity model
    % =========================
    load('Input\Velocal.mat','Vp','topl');
    Vs=Vp/1.73;

    % =========================
    % Fault geometry
    % =========================
    strike=strike1;
    dip=dip1;
    % In my working, I make file by doing for the windows system:
    %       cd C:\MATLAB\MyMechanismFunc\SecondMoment\src\topp
    %       mingw32-make clean
    %       mingw32-make
    % Of course, it's possible that Windows fails to find the runtime library DLLs
    % when launching topp.exe, since it depends on libgfortran-5.dll and libquadmath-0.dll.
    % Thus, we will statically link the MinGW runtime libraries into the topp.exe compilation
    %       LDFLAGS = -static-libgfortran -static-libgcc -static-libquadmath
    [G,takeoffs,taketimes]=Step3_Initial(mlats,mlons,melevs,late,lone,depe,Vp,Vs,topl,phas,strike,dip);
    % where G is inverse operator and takeoffs is take-off angles, taketimes is travel time

    % =========================
    % Semidefinite inversion
    % =========================
    [m2,~,~,Vx,Vy,~,info]=Step3_Inverse(G,d');

    % =========================
    % Derived quantities
    % =========================
    m2=m2(1:6);   % Ditch the dummy variable in the decision vector
    Xvar=[m2(4), m2(5); m2(5), m2(6);]; % Matrix 2X2
    % Misfit results
    Misfit=sum((G*m2-d').^2)./sum(d.^2);
    Misfit_tauc0=sum(abs(2*sqrt(G*m2)-2*sqrt(d')));
    Misfit_tauc1=Misfit_tauc0/sum(2*sqrt(d'));
    % Characteristic quantities
    [U,S,V]=svd(Xvar);
    Lc=2*sqrt(S(1,1)); % or Lc=2*sqrt(max(max(S)));
    Wc=2*sqrt(S(2,2));
    tauc=2*sqrt(m2(1));
    V0=m2(2:3)/m2(1);  % i.e., Vx and Vy
    mV0=sqrt(sum(V0.^2)); % mean Velocity
    Vc=Lc/tauc;
    L0=tauc*mV0;
    ratio=L0/Lc;
    % ratio=mV0/Vc: Directivity ratio, ranges from 0 for a perfectly symmetric bilateral rupture to
    % 1 for a uniform slip unilateral rupture
end
% plots
Fun_PlotEstimate(G,d,m2,mlons,mlats,lone,late,'FigureNumber',14,'MarkerSize',60,...
    'ShowLabels',true);
Fun_PrintEstimate;
fprintf('[3] DONE WITH INVERSION\n');



%%
% >> Sampling uncertainty analysis
% The output may occupy a large portion of the command line window, 
% so please review the above output before executing this module
rng(2026);
isJack=true;
azband=20;
isBoot=true;
Nsample=500;
bconf=0.95;

% Input standardization
G=double(G);
d=double(d(:));
m2=double(m2(:));
if numel(m2)>=7
    m2=m2(1:6);
end
mlons=double(mlons(:));
mlats=double(mlats(:));

% Best-fit derived quantities
best=calDerived(m2,G,d);
% (1) Bootstrap
if isBoot
    [mv0u,mv0l,bound2u,bound2l,Lcu,Lcl,taucu,taucl,boot]=Step4_Bootstrap(...
            G,d,bconf,Nsample,'MakeFigure',true,'FigureNumber',5);
end
% (2) Azimuthal jackknife
if isJack
    jack=Step4_Jackknife(G,d,m2,late,lone,mlats,mlons,azband,'MakeFigure',false,'FigureNumber',6);
end
fclose all;
fprintf('[4] DONE WITH Sampling\n');

