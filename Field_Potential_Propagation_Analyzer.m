% Field Potential Propagation Analyzer
% By Kevin Gray
% Produces electrode activation isochrone plots and calculates and plots 
% conduction velocity x and y components for beat data outputted by the 
% Axion BioSystems Maestro Cardiac Satistics Compiler for a single MEA well. 
% Removes outlier values from conduction velocity and propagation  data. 
% All electrodes must have at least one mutually recorded beat. 

close all; clear all; clc; 

beatTimes = xlsread('MEA_A4.xlsx'); %% read the datafile

%%%%% Format Datafile

for i = 1 : length(beatTimes(:, 1)) %% id and remove NaN values
    if isnan(beatTimes(i, 8)) == 1 
        beatTimes(i, 1) = 0;
        beatTimes(i, 4) = 0;
        beatTimes(i, 8) = 0;
    end
end

first = beatTimes(1, 8); % number of first beat
rawBeats = beatTimes(:, 1);
beatNumbers = beatTimes(:, 8);
beatPeriods = beatTimes(:, 4);
rawBeats = rawBeats(rawBeats~=0);
beatNumbers = beatNumbers(beatNumbers~=0);
beatPeriods = beatPeriods(beatPeriods~=0);
Beats2 = zeros(length(rawBeats), 2);
Beats2(:, 1) = rawBeats;
Beats2(:, 2) = beatNumbers;
beats_2_kill = [];
            
for i = 1 : length(Beats2(:, 2)) - 1; % find non-constituent beats
    if Beats2(i + 1, 2) ~= Beats2(i, 2) + 1 && Beats2(i + 1, 2) ~= first
        beats_2_kill = [(Beats2(i, 2) + 1) beats_2_kill];
    end
end

for i = 1 : length(Beats2(:, 2 - 1)) % remove non-constituent beats
    if ismember(Beats2(i, 2), beats_2_kill) == 1;
        Beats2(i, 1) = 0;
        beatPeriods(i) = 0; 
    end
end

BeatsFinal = Beats2(:, 1);
BeatsFinal = BeatsFinal(BeatsFinal~=0);
beatPeriods = beatPeriods(beatPeriods~=0); 
            

beats = 30 - length(beats_2_kill);

% MEA layout
% 11 .. 18
% :      :
% 81 .. 88
for i = 1:8 % fill in rows
    for j = 1:8 % fill in columns
         MEA(i, j, 1:beats) = BeatsFinal( (j - 1) * beats + (i - 1) * 8 * beats + 1 : j * beats + (i - 1) * 8 * beats);
         MEABeatPeriod(i, j, 1:beats) = beatPeriods( (j - 1) * beats + (i - 1) * 8 * beats + 1 : j * beats + (i - 1) * 8 * beats);
    end
end

MEAzero = zeros(8, 8); 
for beat = 1 : beats
    beatMEA = MEA(:, :, beat);
    MEAzero(:, :) = (MEAzero + (MEA(:, :, beat) - min(beatMEA(:)) * ones(8, 8))) ./ beats;
end

%%%%% Conduction Velocities
[CVx, CVy] = gradient(MEAzero);
CVx = 200 ./ CVx;
CVy = 200 ./ CVy;


% Remove conduction velocities greater than a standard deviation from the
% mean conduction velocity
BeatCVxStd = std2(CVx); 
BeatCVxMean = mean2(CVx); 
BeatCVyStd = std2(CVy); 
BeatCVyMean = mean2(CVy); 
Conduction_Velocity_outliers = 0; 
for j = 1 : 8
    for k = 1 : 8
       if CVx(k, j) > BeatCVxMean + 1.5 * BeatCVxStd || CVx(k, j) < BeatCVxMean - 1.5 * BeatCVxStd
           CVy(k, j) = NaN;
           CVx(k, j) = NaN;
           Conduction_Velocity_outliers = Conduction_Velocity_outliers + 1;  
       elseif CVy(k, j) > BeatCVyMean + 1.5 * BeatCVyStd || CVy(k, j) < BeatCVyMean - 1.5 * BeatCVyStd
           CVx(k, j) = NaN; 
           CVy(k, j) = NaN; 
           Conduction_Velocity_outliers = Conduction_Velocity_outliers + 1;
       end
    end 
end

%%%%% Plot Conduction Velocity Vectors           
[x, y] = meshgrid(8);
quiver(1:8, 1:8, CVx(:, :), CVy(:, :)); %% plot vectors
title('Average Conduction Velocity Vectors')
hold on; 
figure(2)
quiver(1:8, 1:8, CVx(:, :), zeros(8, 8)); %% x component
title('conduction velocity x component')
figure(3)
quiver(1:8, 1:8, zeros(8, 8), CVy(:, :)); %% y component
title('conduction velocity y component')

% Calculate conduction velocity parameters (absolute values of
% conduction velocities used)
Mean_X_Direction_Conduction_Velocity = nanmean(reshape(abs(CVx), 1, numel(CVx)))
Mean_Y_Direction_Conduction_Velocity = nanmean(reshape(abs(CVy), 1, numel(CVy)))
Conduction_Aspect_Ratio = Mean_X_Direction_Conduction_Velocity / Mean_Y_Direction_Conduction_Velocity
Conduction_Velocity_outliers

%%%%% Electrode Activation

% Replace activation times greater than two standard deviations from the
% mean activation time with values from average of x and y linear interpolation
% at that point. 
BeatStd = std2(MEAzero); 
BeatMean = mean2(MEAzero); 
propagation_outliers = 0; 
for j = 1 : 8
    for k = 1 : 8
       if MEAzero(k, j) > BeatMean + 2 * BeatStd || MEAzero(k, j) < BeatMean - 2 * BeatStd
           datarowx = MEAzero(k, :); 
           datarowx(k) = [];
           electrodes = 1 : 8;
           electrodes(k) = []; 
           p = polyfit(electrodes,datarowx, 1);
           datarowy = MEAzero(:, j); 
           datarowy(j) = [];
           electrodes = 1 : 8;
           electrodes(j) = []; 
           p2 = polyfit(electrodes,datarowy', 1); 
           MEAzero(k, j) = (p(1) * k + p(2) + p2(1) * j + p2(2)) / 2; 
           propagation_outliers = propagation_outliers + 1; 
       end
    end 
end

propagation_outliers

figure(4)
[C, h] = contourf(imresize(MEAzero, 200)); 
colormap(jet);
clabel(C,h);
set(gca, 'XTickLabel', [1 2 3 4 5 6 7 8]); 
set(gca, 'YTickLabel', [1 2 3 4 5 6 7 8]);
c = colorbar;
c.Label.String = 'Realative time of activation (milliseconds)';

%%%%% Contour Plot Settings
h.LevelStep = 2 * 10^-5; % Set isochrone spacing
h.LineColor = 'none'; % Comment/uncomment to Turn on/off contour lines
caxis([ 0 .5*10^-3 ]); % Change right value to set time corresponding to red



