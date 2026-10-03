% Build the ray-based design matrix for the second-moment inversion.
% Each row of G corresponds to one station and phase observation. The six
% columns correspond to [mtt, mxt, myt, mxx, mxy, myy] in fault coordinates.
% TOPP is used to calculate the travel time and takeoff angle in a layered
% 1-D velocity model. The input and output files used by TOPP are written to
% MATLAB's current working directory because the Fortran program uses fixed
% file names.
function [G,takeoffAngle,travelTime]=Step3_Initial(stationLat,stationLon,stationElev, ...
    eventLat,eventLon,eventDepth,vpModel,vsModel,layerTop,phase,strike,dip)

if nargin<12
    error('Step3_Initial:NotEnoughInputs','Twelve inputs are required.');
end

stationLat=stationLat(:)';
stationLon=stationLon(:)';
stationElev=stationElev(:)';
phase=upper(char(phase));
phase=phase(:)';
vpModel=vpModel(:)';
vsModel=vsModel(:)';
layerTop=layerTop(:)';

stationCount=numel(stationLat);
if numel(stationLon)~=stationCount || numel(stationElev)~=stationCount || numel(phase)~=stationCount
    error('Step3_Initial:SizeMismatch', ...
        'Station coordinates, elevations, and phase labels must have equal lengths.');
end
if isempty(vpModel) || numel(vpModel)~=numel(vsModel) || numel(vpModel)~=numel(layerTop)
    error('Step3_Initial:VelocityModel', ...
        'Vp, Vs, and layerTop must be nonempty vectors with equal lengths.');
end
if any(~ismember(phase,['P','S']))
    error('Step3_Initial:InvalidPhase','Phase labels must be P or S.');
end
if ~isscalar(eventDepth) || eventDepth<0
    error('Step3_Initial:InvalidDepth','eventDepth must be a nonnegative scalar.');
end

% Locate the platform-specific TOPP executable relative to this file
% According to your own paths
% toppDirectory=fullfile(fileparts(mfilename('fullpath')),'topp');
ToppPath='C:\MATLAB\MyMechanismFunc\SecondMoment\src\topp\';
if ispc
    Toppexe=fullfile(ToppPath,'topp.exe');
else
    Toppexe=fullfile(ToppPath,'topp');
end
% Check exists of topp
if exist(Toppexe,'file')~=2
    error('Step3_Initial:ToppNotFound','TOPP executable was not found: %s. Build it in %s.', ...
        Toppexe,ToppPath);
end
toppCommand=['"',Toppexe,'"'];

% Convert station and event coordinates to a local east-north system.
origin=[mean(stationLat),mean(stationLon)];
stationXY=llh2localxy([stationLat;stationLon;stationElev],origin);
eventXY=llh2localxy([eventLat,eventLon,0]',origin);
stationX=stationXY(:,1);
stationY=stationXY(:,2);
eventX=eventXY(1);
eventY=eventXY(2);

layerCount=numel(vpModel);
G=zeros(stationCount,6);
takeoffAngle=zeros(stationCount,1);
travelTime=zeros(stationCount,1);

for stationIndex=1:stationCount
    stationAzimuth=azimuth(eventLat,eventLon,stationLat(stationIndex),stationLon(stationIndex));
    horizontalDistance=hypot(stationX(stationIndex)-eventX,stationY(stationIndex)-eventY);

    if phase(stationIndex)=='P'
        velocityModel=vpModel;
    else
        velocityModel=vsModel;
    end

    % TOPP reads these fixed names from MATLAB's current working directory.
    distanceInput=horizontalDistance;
    depthInput=eventDepth;
    velocityInput=velocityModel(:);
    topInput=layerTop(:);
    save -ascii toppinputs distanceInput depthInput layerCount velocityInput topInput

    [status,result]=system(toppCommand);
    if status~=0
        error('Step3_Initial:ToppFailed', ...
            'TOPP failed for station %d (status %d): %s',stationIndex,status,strtrim(result));
    end

    outputFile=fopen('toppoutputs','r');
    if outputFile<0
        error('Step3_Initial:ToppOutputMissing', ...
            'TOPP did not create toppoutputs for station %d.',stationIndex);
    end
    travelTime(stationIndex)=fscanf(outputFile,'%g',1);
    takeoffAngle(stationIndex)=fscanf(outputFile,'%g',1);
    fclose(outputFile);

    % Unit ray vector in east, north, down coordinates.
    horizontalRay=sin(takeoffAngle(stationIndex)*pi/180);
    rayVector=[...
        horizontalRay*sin(stationAzimuth*pi/180); ...
        horizontalRay*cos(stationAzimuth*pi/180); ...
        cos(takeoffAngle(stationIndex)*pi/180)];

    % Use the velocity of the layer containing the event depth.
    layerIndex=find(layerTop<=eventDepth,1,'last');
    if isempty(layerIndex)
        layerIndex=1;
    end
    phaseVelocity=velocityModel(min(layerIndex,layerCount));
    slowness=rayVector/phaseVelocity;

    % Rotate from east-north-down to along-strike, normal, down-dip.
    strikeRotation=(strike-90)*pi/180;
    strikeMatrix=[...
        cos(strikeRotation),-sin(strikeRotation),0; ...
        sin(strikeRotation), cos(strikeRotation),0; ...
        0,0,1];
    dipRotation=(dip-90)*pi/180;
    dipMatrix=[...
        1,0,0; ...
        0,cos(dipRotation),sin(dipRotation); ...
        0,-sin(dipRotation),cos(dipRotation)];
    slowness=dipMatrix*strikeMatrix*slowness;

    % G*m predicts the observed temporal second moment at each station.
    G(stationIndex,:)=[...
        1, ...
        -2*slowness(1), ...
        -2*slowness(3), ...
        slowness(1)^2, ...
        slowness(1)*slowness(3), ...
        slowness(3)^2];
end
end
