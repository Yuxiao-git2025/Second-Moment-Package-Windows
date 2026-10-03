% Estimate apparent source-time functions using EGF deconvolution.
%
% This function processes one mainshock-EGF waveform pair for each
% station/component. The user selects the mainshock and EGF windows and
% P- or S-wave arrivals interactively. Projected Landweber deconvolution
% is then used to estimate the apparent source-time function (ASTF).
%
% The returned t2 is the second temporal central moment in s^2:
%
%     t2 = mu^(0,2)
%
% The characteristic duration is:
%
%     tauC = 2*sqrt(t2)
%
% Required external functions:
%     pld.m
%     findt2.m
%     taper.m
%     Fun_defaultAxes.m
%
function [t2,DONE,STF,GFsv,dhatsv,datasv,Tsv,T1sv,epsv,...
    epldsv,tpldsv,t0,t1,PhaseSv]=Step2_Measure(velEGFa,velMSa,niter,npMS,npEGF,dtsva,stasm,compm,pickt2)
% ============ Check input dimensions and basic parameters ============
ns=size(velMSa,1);
Npts=2048;
if size(velEGFa,1)~=ns
    error('velEGFa and velMSa must have the same number of rows.');
end
if numel(npMS)~=ns || numel(npEGF)~=ns || numel(dtsva)~=ns
    error('npMS, npEGF and dtsva must contain one value per station.');
end
if numel(stasm)~=ns || size(compm,1)~=ns
    error('stasm and compm must contain one value per station.');
end
% ============ Select a new measurement set or resume previous work ============
mfile='MEASUREMENTS.mat';
wheremfile='Input\';
% choice=menu('Measure Starts','New','Loading');
choice=menuUI('Measure Starts','New','Loading');

if choice==0
    error('Measurement was cancelled by the user.');
end

if choice==1
    DONE=zeros(ns,1);

    t2=nan(ns,1);
    t1=nan(ns,1);
    t0=nan(ns,1);

    datasv=zeros(ns,Npts);
    GFsv=zeros(ns,Npts);
    STF=zeros(ns,Npts);
    dhatsv=zeros(ns,Npts);

    Tsv=nan(ns,1);
    T1sv=nan(ns,1);
    epsv=nan(ns,1);

    epldsv=nan(ns,Npts);
    tpldsv=nan(ns,Npts);

    PhaseSv=repmat(' ',ns,1);
else
    if exist(mfile,'file')~=2
        error('Cannot find %s.',mfile);
    end
    load(mfile,...
        'DONE','t2','t1','t0','datasv','GFsv','STF','dhatsv',...
        'Tsv','T1sv','epsv','epldsv','tpldsv');
    if exist(mfile,'file')==2
        dataInfo=whos('-file',mfile,'PhaseSv');
        if isempty(dataInfo)
            PhaseSv=repmat(' ',ns,1);
        else
            load([wheremfile,mfile],'PhaseSv');
        end
    else
        PhaseSv=repmat(' ',ns,1);
    end
    if numel(DONE)~=ns
        error('The saved measurement file has a different station number');
    end
end
% Make sure all station-dependent arrays have the expected size
if size(datasv,1)~=ns || size(datasv,2)~=Npts
    error('Saved datasv has an incompatible size');
end

% The input sampling interval is authoritative
dtsv=dtsva(:);
% ============ Process each station/component ============
for i=1:ns
    stationName=strtrim(stasm{i});
    channelName=strtrim(compm(i,:));
    fprintf('[Doing] We now processing: %d-%s-%s\n',i,stationName,channelName);
    mainLength=min(npMS(i),size(velMSa,2));
    egfLength=min(npEGF(i),size(velEGFa,2));
    if mainLength<100 || egfLength<100
        fprintf('[Skip] %s %s: waveform is too short\n',stationName,channelName);
        DONE(i)=0;
        Fun_saveMeasurements();
        continue;
    end

    if ~isfinite(dtsv(i)) || dtsv(i)<=0
        fprintf('Skip %s %s: invalid sampling interval\n',stationName,channelName);
        DONE(i)=0;
        Fun_saveMeasurements();
        continue;
    end
    mainWave=velMSa(i,1:mainLength).';
    egfWave=velEGFa(i,1:egfLength).';

    mainTime=(0:mainLength-1)*dtsv(i);
    egfTime=(0:egfLength-1)*dtsv(i);

    Do=true;
    while Do
        % ============ Display the complete mainshock and EGF waveforms ============
        figure('Name','Mainshock and EGF','NumberTitle','off');
        mscale=max(abs(mainWave));
        egfScale=max(abs(egfWave));
        if mscale>0
            plot(mainTime,mainWave/mscale,'Color',[.2 .2 .2],'LineWidth',0.8);
        else
            plot(mainTime,mainWave,'Color',[.2 .2 .2],'LineWidth',0.8);
        end
        hold on;
        if egfScale>0
            plot(egfTime,egfWave/egfScale,'Color',[0.22 0.37 0.99],'LineWidth',0.8);
        else
            plot(egfTime,egfWave,'Color',[0.22 0.37 0.99],'LineWidth',0.8);
        end
        title([[sprintf('[%d] ',i)],stationName,'.',channelName]);
        xlabel('Time (s)');
        ylabel('Normalized Amp');
        Fun_defaultAxes;

        
        if DONE(i)==1
            stationChoice=menuUI('The Station is already finished, retry?',...
                'Yes','No');
        else
            stationChoice=menuUI('Not finished, go to deconvolution?',...
                'Yes','No');
        end
        close(gcf);


        if stationChoice==0
            error('Measurement was cancelled by the user');
        end
        if stationChoice==2
            Do=false;
            continue;
        end

        % ============ Select the mainshock waveform window ============
        mainZoom=Fun_PickWindow(mainWave,...
            ['[M1] Zoom out for target: ',stationName,' ',channelName]);
        mainInvert=Fun_PickWindow(mainWave(mainZoom(1):mainZoom(2)),...
            ['[M2] Zoom out for inverse: ',stationName,' ',channelName]);
        mainStart=mainZoom(1)+mainInvert(1)-1;
        mainEnd=mainZoom(1)+mainInvert(2)-1;
        mainStart=max(1,min(mainStart,mainLength));
        mainEnd=max(1,min(mainEnd,mainLength));
        if mainEnd<=mainStart
            warning('Invalid mainshock window. Please repeat.'); continue;
        end
        dataRaw=mainWave(mainStart:mainEnd);

        % Limit the data to the maximum allowed length
        dataLength=min(numel(dataRaw),Npts);
        data=zeros(Npts,1);
        data(1:dataLength)=dataRaw(1:dataLength);

        % Pick the P- or S-wave arrival inside the selected mainshock window
        figure('Name','Mainshock Arrival','NumberTitle','off');hold on;
        plot(data(1:dataLength),'Color',[.2 .2 .2],'LineWidth',1.2);
        title('[M3] Select the first arrival for P/S waves: ');
        xlabel('Sample');
        ylabel('Amp');
        Fun_defaultAxes;

        [arrivalSample,~]=ginput(1);
        close(gcf);

        if isempty(arrivalSample) || ~isfinite(arrivalSample)
            warning('No mainshock arrival was selected. Please repeat.');
            continue;
        end

        arrivalSample=round(arrivalSample);
        arrivalSample=max(1,min(arrivalSample,dataLength));

        % Keep a few samples before the picked arrival (e.g., 5~20)
        T1=max(1,arrivalSample-5);

        % Apply a taper only to the selected mainshock data
        data(1:dataLength)=data(1:dataLength).*taper(dataLength,.1);
        datasv(i,:)=data.';


        % ============ Select the EGF arrival ============
        egfZoom=Fun_PickWindow(egfWave,['[E1] Zoom out for inverse: ',stationName,'.',channelName]);

        figure('Name','EGF Arrival','NumberTitle','off');
        egfZoomWave=egfWave(egfZoom(1):egfZoom(2));
        plot(egfZoomWave,'Color',[0.22 0.37 0.99],'LineWidth',1.2);
        title('[E2] Select the first arrival for P/S waves: ');
        xlabel('Sample');
        ylabel('Amp');
        Fun_defaultAxes;

        [egfArrivalLocal,~]=ginput(1);
        close(gcf);

        if isempty(egfArrivalLocal) || ~isfinite(egfArrivalLocal)
            warning('No EGF arrival was selected. Please repeat.');
            continue;
        end
        egfArrivalLocal=round(egfArrivalLocal);
        egfArrivalLocal=max(1,min(egfArrivalLocal,numel(egfZoomWave)));
        egfArrival=egfZoom(1)+egfArrivalLocal-1;

        % ============ Extract a compatible EGF window ============
        % The original implementation uses approximately the data duration
        % after T1 as the EGF duration.
        tarLen=dataLength-T1+1;
        egfEnd=egfArrival+tarLen-1;
        if egfEnd>egfLength
            warning(['The EGF does not contain enough samples after ',...
                'the picked arrival. Please select another window.']);
            continue;
        end
        GF=zeros(Npts,1);
        GF(1:tarLen)=egfWave(egfArrival:egfEnd);
        % Do not taper the EGF onset here. The original method treats the
        % picked EGF arrival as the Green-function start
        GFsv(i,:)=GF.';

        % ============ Projected Landweber deconvolution ============
        try
            [STFs,fitWave,T,fitError,trialTimes,trialErrors]=pld(data,GF,T1,niter);
        catch ME
            warning('PLD failed for %s %s: %s',...
                stationName,channelName,ME.message);
            continue;
        end
        if isempty(STFs)
            warning('PLD returned an empty ASTF. Please repeat.');
            continue;
        end
        STFs=STFs(:);
        fitWave=fitWave(:);
        trialTimes=trialTimes(:);
        trialErrors=trialErrors(:);
        sourceLength=min(numel(STFs),Npts);
        fitLength=min(numel(fitWave),Npts);
        STF(i,:)=0;
        STF(i,1:sourceLength)=STFs(1:sourceLength).';
        dhatsv(i,:)=0;
        dhatsv(i,1:fitLength)=fitWave(1:fitLength).';

        [momentSecond,t1Sample,t0Sample]=findt2(STFs,pickt2);
        t2(i)=momentSecond*dtsv(i)^2; % unit: s^2
        t1(i)=t1Sample;
        t0(i)=t0Sample;

        Tsv(i)=(T-T1)*dtsv(i);
        T1sv(i)=T1*dtsv(i);
        epsv(i)=fitError;

        epldsv(i,:)=nan;
        tpldsv(i,:)=nan;

        trialLength=min(numel(trialErrors),Npts);
        epldsv(i,1:trialLength)=trialErrors(1:trialLength).';

        trialLength=min(numel(trialTimes),Npts);
        tpldsv(i,1:trialLength)=trialTimes(1:trialLength).';




        % ============ Show the result of measurement summary ============
        figure('Name','Measurement summary','NumberTitle','off');
        tiledlayout(2,2,"TileSpacing","compact","Padding","compact");
        nexttile;% Plot-1
        dataScale=max(abs(data));
        gfScale=max(abs(GF));
        if dataScale>0
            plot((0:Npts-1)*dtsv(i),data/dataScale,'Color',[.2 .2 .2],'LineWidth',1.2);
        else
            plot((0:Npts-1)*dtsv(i),data,'Color',[.2 .2 .2],'LineWidth',1.2);
        end
        hold on;
        if gfScale>0
            plot((0:Npts-1)*dtsv(i)+T1*dtsv(i),GF/gfScale,'Color',[0.22 0.37 0.99],'LineWidth',1.2);
        else
            plot((0:Npts-1)*dtsv(i)+T1*dtsv(i),GF,'Color',[0.22 0.37 0.99],'LineWidth',1.2);
        end
        xlabel('Time (s)');
        ylabel('Normalized Amp');
        Fun_defaultAxes;
        lastIdx=find(data~=0,1,'last');
        if isempty(lastIdx)
            lastIdx=1;
        end
        tEnd=(lastIdx-1)*dtsv(i);
        xlim([0, tEnd]); % Reset the x-axis range
        title('Origin Waveforms');



        nexttile; % Plot-2
        sourceScale=max(abs(STF(i,:)));
        if ~isfinite(sourceScale) || sourceScale<=0
            sourceScale=1;
        end
        tAxis=(0:Npts-1)*dtsv(i);
        stfPlot=STF(i,:)/sourceScale;
        patch('XData',[tAxis,fliplr(tAxis)],...
            'YData',[stfPlot,zeros(size(stfPlot))],...
            'FaceColor',[0.00 0.90 0.90],...
            'FaceAlpha',0.2,...
            'EdgeColor','none');
        hold on;
        plot(tAxis,stfPlot,'Color',[0.00 0.90 0.90],'LineWidth',1.8);
        if isscalar(t0Sample) && isfinite(t0Sample)
            t0Index=round(t0Sample);

            if t0Index>=1 && t0Index<=size(STF,2)
                t0Time=(t0Sample-1)*dtsv(i);
                t0Value=STF(i,t0Index)/sourceScale;

                plot(t0Time,t0Value,...
                    'Marker','pentagram',...
                    'MarkerFaceColor',[0.88 0.23 0.23],...
                    'MarkerEdgeColor',[0.88 0.23 0.23],...
                    'MarkerSize',20);
            end
            xMax=3*max(t0Sample-1,1)*dtsv(i);
        else
            xMax=max(tAxis);
        end
        title(['ASTF (\tau=',...
            num2str(2*sqrt(max(t2(i),0)),4),' s)']);
        xlabel('Time (s)');
        ylabel('Normalized amplitude');
        xlim([0,max(xMax,dtsv(i))]);
        Fun_defaultAxes;


        nexttile;% Plot-3
        plot(trialTimes*dtsv(i),trialErrors,'Color',[.5 .5 .5],'LineWidth',2.4);
        hold on;
        if isfinite(T)
            [~,bestIndex]=min(abs(trialTimes-T));
            if ~isempty(bestIndex)
                plot(trialTimes(bestIndex)*dtsv(i),trialErrors(bestIndex),'Marker','pentagram', ...
                    'MarkerFaceColor',[0.88 0.23 0.23],'MarkerEdgeColor',[0.88 0.23 0.23],'MarkerSize',20);
            end
        end
        title('PLD Misfit');
        xlabel('ASTF duration (s)');
        ylabel('Misfit');
%         xticks(0:1:20)
        Fun_defaultAxes;



        nexttile;% Plot-4
        plot((0:Npts-1)*dtsv(i),data,'Color',[.2 .2 .2],'LineWidth',1.2);
        hold on;
        plot((0:Npts-1)*dtsv(i),dhatsv(i,:),'Color',[1 .1 .1],'LineWidth',1.2);
        title('Seismogram Fitting');
        xlabel('Time (s)');
        ylabel('Amp');
        Fun_defaultAxes;
        xlim([0, tEnd]);

        set(gcf,"Position",[800, 300, 1100, 700])

        % ============ Decide whether to save the result ============
        saveChoice=menuUI('Save it as',...
            'P-wave',...
            'S-wave',...
            'Retry',...% Restart
            'Abort');% cancel this channel
        close(gcf);
        if saveChoice==0
            error('Measurement was cancelled by the user.');
        elseif saveChoice==1
            DONE(i)=1;
            PhaseSv(i)='P';
            Do=false;
        elseif saveChoice==2
            DONE(i)=1;
            PhaseSv(i)='S';
            Do=false;
        elseif saveChoice==3
            DONE(i)=0;
            Do=true;
        elseif saveChoice==4
            DONE(i)=0;
            Do=false;
        end
    end
    fprintf('[Save] to %s%s\n\n',wheremfile,mfile);
    Fun_saveMeasurements();
end
fprintf('[Finished] Make Measurements of all stations \n\n');


    % Local functions
    function selectedWindow=Fun_PickWindow(waveform,windowTitle)
        waveform=waveform(:);

        figure('Name','Select waveform window','NumberTitle','off');
        plot(waveform,'Color',[.2 .2 .2],'LineWidth',1.2);
        title(windowTitle);
        xlabel('Sample');
        ylabel('Amp');
        Fun_defaultAxes;

        [x1,~]=ginput(1);
        xline(x1,'LineWidth',1.6,'Color',[0.65 0.65 0.65],'LineStyle',':');
        [x2,~]=ginput(1);
        xline(x2,'LineWidth',1.6,'Color',[0.65 0.65 0.65],'LineStyle',':');
        close(gcf);

        if isempty(x1) || isempty(x2) || ...
                ~isfinite(x1) || ~isfinite(x2)
            error('Waveform window selection was cancelled');
        end
        selectedWindow=round(sort([x1,x2]));
        selectedWindow(1)=max(1,min(selectedWindow(1),numel(waveform)));
        selectedWindow(2)=max(1,min(selectedWindow(2),numel(waveform)));
        if selectedWindow(2)<=selectedWindow(1)
            error('The selected waveform window is invalid.');
        end
    end
    
    function Fun_saveMeasurements()
    
    save([wheremfile,mfile],...
        'DONE','PhaseSv','Npts','t2','t1','t0','datasv','dtsv',...
        'dhatsv','GFsv','STF','Tsv','T1sv','epsv',...
        'epldsv','tpldsv','-v7');
    
    % Update the MATLAB base workspace for real-time inspection.
    assignin('base','DONE',DONE);
    assignin('base','PhaseSv',PhaseSv);
    assignin('base','t2',t2);
    assignin('base','t1',t1);
    assignin('base','t0',t0);
    assignin('base','datasv',datasv);
    assignin('base','dtsv',dtsv);
    assignin('base','GFsv',GFsv);
    assignin('base','STF',STF);
    assignin('base','dhatsv',dhatsv);
    assignin('base','Tsv',Tsv);
    assignin('base','T1sv',T1sv);
    assignin('base','epsv',epsv);
    assignin('base','epldsv',epldsv);
    assignin('base','tpldsv',tpldsv);
    
    end



end


% function choice=menuUI(titleStr, opt1, opt2)
%     fig = uifigure('Name', 'Menu', ...
%                    'Position', [400 900 550 250], ...
%                    'Resize', 'off','Color',[1 1 1]);
%     uilabel(fig, 'Text', titleStr, ...
%             'FontSize', 20, 'FontWeight', 'bold', ...
%             'Position', [30 130 300 40], ...
%             'HorizontalAlignment', 'center');
%     choice = 0;
%     uibutton(fig, 'Text', opt1, 'FontSize', 18, ...
%              'Position', [50 40 110 50], ...
%              'ButtonPushedFcn', @(src,evt) setChoiceAndClose(1));
%     uibutton(fig, 'Text', opt2, 'FontSize', 18, ...
%              'Position', [200 40 110 50], ...
%              'ButtonPushedFcn', @(src,evt) setChoiceAndClose(2));
%     uiwait(fig);
%     function setChoiceAndClose(val)
%         choice = val;
%         delete(fig);
%     end
% end
function choice = menuUI(titleStr, varargin)
%   choice = menuUI(titleStr, opt1, opt2, ..., optN)
%   Return value: 1..N corresponds to each button; 0 means the user closed the window directly
    nOpt = numel(varargin);
    if nOpt == 0
        error('menuUI: At least one option is required.');
    end
    % ===== Compute window size dynamically based on number of buttons =====
    btnW   = 150;      % Button width
    btnH   = 55;       % Button height
    gapX   = 50;       % Horizontal gap between buttons
    gapY   = 25;       % Vertical gap between buttons
    margin = 30;       % Margin around the edges
    maxPerRow = 2;                          % Max number of buttons per row
    nRow = ceil(nOpt / maxPerRow);          % Total number of rows
    nCol = min(nOpt, maxPerRow);            % Max number of columns per row
    % Content area width
    contentW = nCol*btnW + (nCol-1)*gapX;
    figW = contentW + 2*margin;
    % Height = title + button area + top/bottom margins
    titleH = 60;
    figH = margin + nRow*btnH + (nRow-1)*gapY + titleH + margin;
    % ===== Create figure =====
    fig = figure('Name', 'Menu', ...
                 'NumberTitle', 'off', ...
                 'MenuBar', 'none', ...
                 'ToolBar', 'none', ...
                 'Resize', 'off', ...
                 'Color', [1 1 1], ...
                 'Position', [400, 600, figW, figH], ...
                 'WindowStyle', 'normal'); % The WindowStyle should be normal rather than modal

    % ===== Title (horizontally centered, near the top) =====
    uicontrol(fig, 'Style', 'text', ...
              'String', titleStr, ...
              'FontSize', 14, 'FontWeight', 'bold', ...
              'BackgroundColor', [1 1 1], ...
              'Position', [margin, figH - margin - titleH, figW - 2*margin, titleH]);

    choice = 0;

    % ===== Create buttons one by one =====
    for k = 1:nOpt
        % Compute the row and column of this button
        row = ceil(k / maxPerRow);           % Which row (from top to bottom)
        col = mod(k-1, maxPerRow) + 1;       % Which column in this row

        % Actual number of buttons in this row (the last row may be partially filled)
        nColThisRow = min(maxPerRow, nOpt - (row-1)*maxPerRow);
        rowContentW = nColThisRow*btnW + (nColThisRow-1)*gapX;
        rowStartX   = (figW - rowContentW)/2;   % Starting X for horizontal centering of this row

        btnX = rowStartX + (col-1)*(btnW + gapX);

        % bottom is measured from the bottom up: the last row is at the bottom
        btnY = margin + (nRow - row)*btnH + (nRow - row)*gapY;

        uicontrol(fig, 'Style', 'pushbutton', ...
                  'String', varargin{k}, ...
                  'FontSize', 18, ...
                  'Position', [btnX, btnY, btnW, btnH], ...
                  'Callback', @(src, evt) setChoiceAndClose(k));
    end

    uiwait(fig);

    % ===== Nested callback =====
    function setChoiceAndClose(val)
        choice = val;
        delete(fig);
    end
end