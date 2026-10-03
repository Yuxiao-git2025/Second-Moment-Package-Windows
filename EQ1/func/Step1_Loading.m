% Build the SJFsetup input without Antelope.
% XA,Yu 2026/09/29
% This script reproduces the variables written by runSJFsetup.m, but uses:
%   1) the flat Antelope text tables SJFexdb.site and SJFexdb.origin only
%      as ordinary text metadata files;
%   2) the bundled MATLAB MiniSEED reader rdmseed.m;
%   3) the local MSdata and EGFdata directories.
%
clearvars;
clc;
%% Paths and output
PATH = 'C:\MATLAB\MyMechanismFunc\SecondMoment\SJFex\';

outputFile = fullfile(PATH,'Input\Data.mat');
msDir = fullfile(PATH,'MSdata');
egfDir = fullfile(PATH,'EGFdata');
siteFile = fullfile(PATH,'Input\SJFexdb.site');
originFile = fullfile(PATH,'Input\SJFexdb.origin');


assert(exist(siteFile,'file')==2, 'Missing site metadata file: %s', siteFile);
assert(exist(originFile,'file')==2, 'Missing origin metadata file: %s', originFile);
assert(exist(msDir,'dir')==7, 'Missing mainshock directory: %s', msDir);
assert(exist(egfDir,'dir')==7, 'Missing EGF directory: %s', egfDir);

%% Event and processing parameters

MSevid = 1621;  % MS ID
EGFevid = 1755; % EGF ID
emag = 5.1;     % Magnitude
dointe = 0;     % Integral control
dofilt = 1;     % Filter control
fmin = 0.25;    % Filter bandwidth
fmax = 15;
strike1 = 307;  % Nodal plane solution
dip1 = 83;
strike2 = 216;
dip2 = 82;

% Station/component selection
stasm = {'TRAN'; 'BALD'; 'IWR'; 'LVA2'; 'SND'; 'HSSP'; 'CRY'; ...
         'TRO'; 'SETM'; 'BZN'; 'BSAP'; 'RRSP'; 'MGD'; 'MGE'; 'SLR'; ...
         'SAL'; 'DEV'; 'PSD'; 'MSC'; 'SNO'; 'PFO'};
compm = ['HNN'; 'EHN'; 'HNN'; 'HNE'; 'HNE'; 'HNE'; 'HNN'; ...
         'HNE'; 'HNN'; 'HHN'; 'HNE'; 'HHN'; 'HNE'; 'HNE'; 'HHN'; ...
         'HHN'; 'HHZ'; 'HHZ'; 'HHZ'; 'HHZ'; 'HNZ'];
% Notice stasm and compm must have the same number of rows
ns = numel(stasm);
nc = numel(compm);
fprintf('Have selected %d components of waves\n',ns)

%% Read station and event metadata
[siteNames, siteLat, siteLon] = readSiteTable(siteFile);
[late, lone, depe] = readOriginTable(originFile, MSevid);
% Confirm the EGF event exists, although its location is not required below.
readOriginTable(originFile, EGFevid);

slat = nan(ns,1);
slon = nan(ns,1);
waitHandle=waitbar(0,'Processing');
for i = 1:ns
    k = find(strcmp(siteNames, strtrim(stasm{i})), 1, 'last');
    if isempty(k)
        error('Station %s is absent from %s.', strtrim(stasm{i}), siteFile);
    end
    slat(i) = siteLat(k);
    slon(i) = siteLon(k);
    waitbar(i/ns,waitHandle,sprintf('[Now] %.0f',i/ns*100));
    pause(0.02);
end
pause(1);
close(waitHandle);

%% Read and process waveform pairs
% The original script pads each station row to the maximum waveform length.
npMS = zeros(ns,1);
npEGF = zeros(ns,1);
dtsv = nan(ns,1);
msData = cell(ns,1);
egfData = cell(ns,1);

for i = 1:ns
    station = strtrim(stasm{i});
    component = strtrim(compm(i,:));

    msFile = findWaveform(msDir, station, component);
    if isempty(msFile)
        error('[%d] No MS.MiniSEED file found for %s %s in %s',i, station, component, msDir);
    end
    egfFile = findWaveform(egfDir, station, component);
    if isempty(egfFile)
        warning('[%d] No EGF.MiniSEED file found for %s %s in %s',i, station, component, egfDir);
    end

    % The supplied setup integrates rows whose second component letter is N (N is acceleration).
    % The channel code consists of three uppercase letters or digits, where the first digit 
    % represents the frequency band code, the second digit represents the
    % instrument type code, and the third digit represents the channel direction code
    if component(2)=='H'
        integrateMS = false;
    elseif component(2)=='N'
        integrateMS = true;
    else
        integrateMS = false;
    end

    [msTrace, dt] = readTrace(msFile, dofilt, fmin, fmax, integrateMS);
    if isempty(egfFile)
        egfTrace = zeros(1,0);
    else
        [egfTrace, dtEGF] = readTrace(egfFile, dofilt, fmin, fmax, integrateMS);
        if abs(dt-dtEGF) > max(1e-10, 1e-6*dt)
            error('Sampling interval mismatch at %s %s: MS=%g, EGF=%g.', station, component, dt, dtEGF);
        end
    end

    dtsv(i) = dt;
    msData{i} = msTrace(:).';
    egfData{i} = egfTrace(:).';
    npMS(i) = numel(msData{i});
    npEGF(i) = numel(egfData{i});
end
% Extend the dimension by largest scalar dimension
maxMS = max(npMS);
maxEGF = max(npEGF);
velMS = zeros(ns,maxMS);
velEGF = zeros(ns,maxEGF);
for i = 1:ns
    velMS(i,1:npMS(i)) = msData{i};
    velEGF(i,1:npEGF(i)) = egfData{i};
end
fprintf('Initial process is completed \n');

%% Save a drop-in replacement for the Antelope-generated file
% Total 17 variables required
save(outputFile, 'compm', 'depe', 'dip1', 'dip2', 'dtsv', 'late', 'lone', ...
    'npEGF', 'npMS', 'slat', 'slon', 'stasm', 'strike1', 'strike2', ...
    'velEGF', 'velMS', '-v7');

fprintf('Wrote %s\n', outputFile);
fprintf('Stations: %d | MS samples: %d~%d | EGF samples: %d~%d\n', ...
    ns, min(npMS), max(npMS), min(npEGF), max(npEGF));
fprintf('Mainshock: lat %.5f, lon %.5f, depth %.3f km\n', late, lone, depe);
%%








% Local helpers
function filePath = findWaveform(folder, station, component)
    patt = sprintf('%s.*.%s.mseed', station, component);
    files = dir(fullfile(folder,patt));
    if isempty(files)
        % A permissive fallback supports filenames with a nonstandard network code.
        files = dir(fullfile(folder,sprintf('%s*%s*.mseed',station,component)));
    end
    if isempty(files)
        filePath = '';
    elseif numel(files)==1
        filePath = fullfile(folder,files(1).name);
    else
        names = {files.name};
        hasStation = ~cellfun('isempty',strfind(names,[station,'.']));
        hasComponent = ~cellfun('isempty',strfind(names,['.',component,'.']));
        exact = find(hasStation & hasComponent,1);
        if isempty(exact)
            error('Multiple candidate MiniSEED files for %s %s.', station, component);
        end
        filePath = fullfile(folder,files(exact).name);
    end
end

function [trace, dt] = readTrace(filePath, dofilt, fmin, fmax, integrateTrace)
    X = rdmseed(filePath);
    if isempty(X)
        error('No records decoded from %s.', filePath);
    end
    t = cat(1,X.t); %#ok<NASGU>
    d = cat(1,X.d);
    if isempty(d)
        error('No samples decoded from %s.', filePath);
    end
    sampr = X(1).SampleRate;
    dt = 1./sampr;
    if dofilt
        if fmax >= sampr/2
            error('fmax=%g Hz is at or above Nyquist=%g Hz for %s.', fmax, sampr/2, filePath);
        end
        [B,A] = butter(4,2*[fmin fmax]/sampr);
        d = d-mean(d);
        d = detrend(d);
        d = d.*taper(length(d),.01);
        d = filtfilt(B,A,d);
    else
        d = d-mean(d);
    end
    if integrateTrace
        trace = inte(d-mean(d),dt);
    else
        trace = d-mean(d);
    end
    trace = trace(:).';
end

function [names, lat, lon] = readSiteTable(filePath)
    fid = fopen(filePath,'r');
    assert(fid>0,'Cannot open %s.',filePath);
    names = {}; lat = []; lon = [];
    cleaner = onCleanup(@() fclose(fid));
    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        tok = regexp(line,'^\s*(\S+)\s+\S+\s+\S+\s+([-+]?\d+(?:\.\d*)?)\s+([-+]?\d+(?:\.\d*)?)','tokens','once');
        if isempty(tok), continue; end
        names{end+1,1} = tok{1}; %#ok<AGROW>
        lat(end+1,1) = str2double(tok{2}); %#ok<AGROW>
        lon(end+1,1) = str2double(tok{3}); %#ok<AGROW>
    end
end

function [lat, lon, depth] = readOriginTable(filePath, evidWanted)
    fid = fopen(filePath,'r');
    assert(fid>0,'Cannot open %s.',filePath);
    lat=[]; lon=[]; depth=[];
    cleaner = onCleanup(@() fclose(fid));
    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        tok = regexp(line,['^\s*([-+]?\d+(?:\.\d*)?)\s+([-+]?\d+(?:\.\d*)?)\s+([-+]?\d+(?:\.\d*)?)\s+' ...
            '[-+]?\d+(?:\.\d*)?\s+\S+\s+(\d+)'],'tokens','once');
        if isempty(tok), continue; end
        if str2double(tok{4})==evidWanted
            lat=str2double(tok{1}); lon=str2double(tok{2}); depth=str2double(tok{3});
            return;
        end
    end
    error('Event evid %d is absent from %s.', evidWanted, filePath);


end

