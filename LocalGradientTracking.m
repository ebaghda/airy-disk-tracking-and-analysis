%% Local Gradient Detection for Particle Tracking
% Select data THIS IS CURRENT VERSION AS OF 09.13.2024
clc; clear; close all; startTic = tic; 
cd("C:\Users\baghd\OneDrive\Desktop\msd_local_gradient"); %%%%%%% SELECT LOCAL FOLDER CONTAINING THIS CODE %%%%%%%
dateFile = "I:\2024.08.29_passive"; %select the experiment date you wish to track with folder structure date>experiement>videos
generateFigures = true; %toggle figure generation
saveData = false; %toggle results file generation
exportTiffs = false; %toggle .tiff figure generation
expDirs = listexpdirs(dateFile); 
head(expDirs)
allMSDData = cell([]); %preallocate
for eNum =1:1%numel(expDirs) %%%%%%%%%%%%%%% loop through experiment folders %%%%%%%%%%%%%%%%
    [nd2Files, nd2Count, fileFullPaths] = listnd2s(expDirs,eNum); 
    msdData = cell([]); %preallocate
    for vidIdx = 1:1%height(nd2Files) %%%%%%%%%%%%%%%%%% loop over nd2 files %%%%%%%%%%%%%%%
        bfrObj = BioformatsImage(char(fileFullPaths(nd2Count(vidIdx))));
        maxFrames = round(bfrObj.sizeT);
        wb = waitbar(0, "Importing ND2 file");
        for i=1:maxFrames %load the first set of frames into matlab
            wb = waitbar(i/maxFrames, wb, "Importing ND2 file");
            try
                scaleFactor=1;% bins video data. Default is 1
                imStack(:,:,i) = double(imresize(getPlane(bfrObj, 1, 1, i), scaleFactor));
            catch ME
                warning("File read error encountered")
                disp(ME)
                imStack=zeros([1,1,1]);
                return
            end
        end
        close(wb);
        if generateFigures
            figure(); imshowpair(imStack(:,:,1),imStack(:,:,maxFrames))
        end

        % Perform object detection using local gradient algorithm
        R=5; % set window size default=3
        thrsize=3; % threshold value default=3
        thrtype= 'topfraction'; % threshold type default=topfraction alt = topvalue
        epsilon=20; % neighborhood search radius (see DBSCAN) default 3
        minpts=4; %  minimum number of neighbors to identify a core point (see DBSCAN)
        %set trajectory linking parameters
        dc=4; % max distance from detected particle to look for linked particles in other frames
        dfr=2; % number of frames to look for linked particles default=3
        Nfr_min=100; % minimum number of frames in trajectory default = 10
        %imStack = zeros(size(img));
        %imStack(:,:,:)=double(img(:,:,:)); % convert to preferred format
        wb = waitbar(0, "Detecting particles");
        LGradTime = tic; %start localization
        coord{maxFrames,1}=[];
        for i=1:maxFrames
            wb = waitbar(i/maxFrames, wb, "Detecting Local Gradients");
            Im=imStack(:,:,i);
            if range(imStack(:,:, i), "all") > 500 % detects overexposed frames
                xyarray = LocalGradient.local_gradient_multi(Im,R,thrsize,thrtype,epsilon,minpts);
                coord{i,1} = xyarray;
                % if length(coord{1,1}) > 40 % if more than 20 particles are detected in one frame
                %     warning("Too many particles detected")
                %     break %end loop
                % end
            else
                i = i+75;
            end
        end
        close(wb)
        disp(['Average execution time per image: ', num2str(toc(LGradTime)/maxFrames*1000,'%.1f'), ' ms' ])

        % Link positions into trajectories
        for i1 = 1:height(coord)
            if isempty(coord{i1,1})  %this fixes a bug
                coord{i1,1} = [0.0001, 0.0001]; %this fixes a bug
            end
        end
        T_fr = LocalGradient.detect_trj2D_eb(coord,dc,dfr,Nfr_min);


        % Show an image with detected particles and trajectories
        if generateFigures
            Im_N = maxFrames; % frame number to show
            f=figure(); imshow(imStack(:,:,Im_N),[]);colormap(parula);hold on;
            plot(coord{Im_N}(:,1),coord{Im_N}(:,2),'ro', 'MarkerSize',20), colorbar % plot circles around particles
            for i=1:height(T_fr)
                plot(T_fr.xy{i}(:,1),T_fr.xy{i}(:,2),'-',"Color","#7E2F8E", "LineWidth", 1.5) % add trajectories
            end
            hold off
            if exportTiffs
                exportgraphics(f, sprintf("trajectories00%i.tiff", vidIdx), "Resolution", 1200);
            end
        end

        tracks = cell(height(T_fr),1);
        for i=1:height(T_fr) %store the trajectories in a tracks cell array
            tracks(i) = {[T_fr.frames{i}(:)./100 T_fr.xy{i}(:,1)*bfrObj.pxSize(1) T_fr.xy{i}(:,2)*bfrObj.pxSize(1)]};
        end

        if generateFigures
            [f, ax] = newfig(); hold on %plot the tracks
            for i=1:height(T_fr) %store the trajectories in a tracks cell array and plot
                tracks(i) = {[T_fr.frames{i}(:)./100 T_fr.xy{i}(:,1)*bfrObj.pxSize(1)/scaleFactor T_fr.xy{i}(:,2)*bfrObj.pxSize(1)/scaleFactor]}; %store trajectories in tracks cell array for msdanalyzer
                if numel(T_fr.xy{i,1}(:,1))>10
                    scatter(T_fr.xy{i,1}(:,1), T_fr.xy{i,1}(:,2), ".")
                end
            end
            xlabel("x-position (px)"); ylabel("y-positon (px)")
            daspect([1 1 1])
            ax.XLim=[0 bfrObj.width]; ax.YLim=[0 bfrObj.height];
            ax.Box = "on"; hold off
        end
        if exportTiffs
            exportgraphics(f, sprintf("tracks00%i.tiff", vidIdx), "Resolution", 1200);
        end

        clearvars ma
        ma = msdanalyzer(2, 'µm', 's');
        try
            ma = ma.addAll(tracks);
            ma = ma.computeMSD;
        catch MA
            warning("Issue in adding tracks. Skipping video")
            disp(MA)
            continue
        end

        % f = figure(); %plot MSD curves
        % ax = axes();
        % ax.Title.String="MSD Curves with 95% CI";
        % ax.set("FontSize",12);
        % ax.Box="on";
        % ax.XGrid="on";
        % ax.YGrid="on";
        % ax.XLabel.String = "Lag Time (s)";
        % ax.YLabel.String = "Mean-Squared Displacement (µm²)";
        % hold on
        clear fitV
        fitV=[];
        for i= 1:height(tracks)
            if isobject(ma) && ma.msd{i}(1,4)> 5 && var(ma.tracks{i}(2,:)) > 0.3
                %%errorbar(ma.msd{i}(:,1),ma.msd{i}(:, 3)); %STD
                %errorbar(ma.msd{i}(1:11,1),ma.msd{i}(1:11, 2), ma.msd{i}(1:11, 3)./sqrt(ma.msd{i}(1:11, 4)).*tinv([0.05], ma.msd{i}(1:11, 4)),"k.", "MarkerSize",0.001); %95% CI
                [fitObj, gof] = FitMSD(ma.msd{i}(1:11, 1),ma.msd{i}(1:11, 2));
                fitV = [fitV, coeffvalues(fitObj)];
                %plot(fitObj, "g-");
                %scatter(ma.msd{i}(1:11,1),ma.msd{i}(1:11, 2), "filled")
                %legend off
            end
        end
        %hold off

        msdData(1, vidIdx) = {nd2Files.name(vidIdx)};
        msdData(2,vidIdx)={ma};
        msdData(3, vidIdx) = {"fit velocity"};
        msdData(4, vidIdx) = {fitV};
        fprintf("Experiment %1$i video %2$i contains %3$i track(s) with the following velocities:", eNum, vidIdx, height(tracks));
        disp(fitV)
        clearvars -except allMSDData fileFullPaths nd2Files nd2Count vidIdx eNum msdData isfolder dateFile dateDir expDirs startTic generateFigures saveData exportTiffs
    end % end looping over videos

    allMSDData(eNum) = {msdData};
    if saveData
        cd(dateFile)
        %ProcessedDateName = sprintf("MSDsInCellsAll_240428_E%d.mat", eNum)
        %ProcessedDateName = sprintf("MSDsInCellsAll_brightfield_240430_E%d.mat", eNum)

        save(ProcessedDateName)
        fprintf("Saving data for experiment number %d", eNum)
        cd("C:\Users\baghd\Desktop\msd_local_gradient")
    end
end %end looping over experiment folders
toc(startTic)
if ~exist("MException","class")
    disp("Code completed without errors!")
end

velocity = collectVs(allMSDData)

open velocity

vIndexed = []; vNameIndexed = [];
for i1 = 1:numel(velocity)
    vCurrent = zeros(height(velocity{i1}), 1);
    vCurrent(:,1) = velocity{i1};
    vNameCurrent = string(zeros(height(velocity{i1}), 1));
    vNameCurrent(:,1) = repmat(expDirs(i1), [height(vCurrent),1]);
    vIndexed  = [vIndexed; vCurrent];
    vNameIndexed = [vNameIndexed; vNameCurrent];
end
open vIndexed
open vNameIndexed


vMatrix  = zeros(300);
for i1 = 1:numel(velocity)
    filteredV = velocity{i1};
    filteredV = filteredV(filteredV>0.1);
    vMatrix(1:numel(filteredV),i1) = filteredV;
end
open vMatrix

%V = readtable("L:\2024.03.05\driftVelocities.csv")

% [f, ax] = newfig(); hold on
% title("Increasing DMSO decreases velocity");
% scatter(ax, V.mf_DMSO, V.Drift_Velocity, 40, V.mf_DMSO, 'filled', 'o', 'MarkerEdgeColor','k');
% xlabel("Mole Fraction DMSO (-)"); ylabel("Drift Velocity (µm/s)");
% xlim([-0.02 0.3]);
% hold off
% %exportgraphics(f, "L:\2024.03.05\scatterVinDMSO.tiff", "Resolution",1200);
% [f, ax] = newfig(); hold on
% ax.PlotBoxAspectRatio = [5 3 1];
% ax.LineWidth = 1.25;
% title("Increasing DMSO decreases velocity");
% boxplot(ax, V.Drift_Velocity,categorical(V.mf_DMSO), "BoxStyle","outline", "Notch","marker", ...
%     "Symbol", 'o', "MedianStyle","line", "OutlierSize",4, "FactorDirection","list", "Widths",0.7, ...
%     "FullFactors","on", "ExtremeMode","clip");
% xlabel("Mole Fraction DMSO"); ylabel("Drift Velocity (µm/s)");
% hold off
% %exportgraphics(f, "L:\2024.03.05\boxVinDMSO.tiff", "Resolution",1200);

function velocity = collectVs(allMSDData)
velocity = {};
for i1 = 1:width(allMSDData)
    V = [];
    for i2 = 1:width(allMSDData{i1})
        V = [V, allMSDData{i1}{4, i2}];
    end
    velocity{i1} = V';
end

end
