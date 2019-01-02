% plotGEFSfreezingLevels.m
% downloads and plots plan view maps of GEFS probabilistic freezing level
% heights relative to California terrain at different forecast hours 
%
% Brian Henn, CW3E/SIO/UCSD, February 2017
% 
% v2 created in May 2017 to include QPE in plots, extend to HUC-8s, show
% basin-average values on main map and detailed info on insets
% v3 created in November 2017 to create plots across the Western US
% compatible with the CW3E AR decision support website for DWR
% v4 created in January 2017 to use WPC precipitation
% _skyriver version designed to run on CW3E server

clear; close all; tic;

%% inputs

% processed freezing level data file from NOMADS
freezingLevelFile = '../out.nc';
precipPath = '../';

% North America 1km DEM file
DEMfile = './MATfiles/NorthAmerica1kmDEM.mat';

% HUC8 basin shapefile
% westCoastHUC8Shapefile = '../shapefiles/WestCoastHUC8s_clipped.shp';

% state boundary shapefile
% stateShapes = '../shapefiles/US_Can_Mex_provinces.shp';

% basins and masks
maskFile = './MATfiles/WestCoastHUC8Masks.mat';

% plotting bounding box
boundingBox = [-128 -114 32 51];

% colormap file
cmapFile = './MATfiles/RainSnowColorMap.mat';

% divergent colormap
cmap3file = 'MATfiles/DivergentCMap.mat';

% forecast initialization time
initTime = datetime(2018,12,22,12,0,0);
% initTime = dateshift(datetime(datestr(now)),'start','hour');

% emsemble dimensions and timesteps
nEns = 21; % ensemble members
inct = 6; % forecast steps to plot (hours)
lenf = 192; % length of GEFS forecast period (hours)
incf = 3; % increment of GEFS forecast (hours)
incp = 6; % increment of WPC precipitation (hours)
lenp = 168; % length of WPC precipitation forecast (hours)

% figure output directory
outputDir = '../outputFiles';

% figure archive directory
archiveDir = '/data/downloaded/Forecasts/ARPortal_Archive/gefs/RainSnowTool';

% offset for Z0C to rain vs. snow (m)
zOffset = 180;

%% create time array

nt = round(lenf/inct) + 1;
time = initTime:hours(inct):(initTime + hours(lenf));
ntP = round(lenp/incp);

%% subset and load freezing level forecasts

lat = ncread(freezingLevelFile,'lat');
lon = ncread(freezingLevelFile,'lon');
lon(lon > 180) = lon(lon > 180) - 360;
goodLon = find(lon >= boundingBox(1) & lon <= boundingBox(2));
goodLat = find(lat >= boundingBox(3) & lat <= boundingBox(4));
[LonGEFS,LatGEFS] = meshgrid(lon(goodLon),lat(goodLat));
Z0C = NaN(nEns,nt,size(LatGEFS,1),size(LatGEFS,2));
for i = 1:nEns
    for j = 1:nt
        Z0C(i,j,:,:) = double(ncread(freezingLevelFile,'Z0C',...
            [goodLon(1) goodLat(1) j*round((inct/incf)) - 1 i],[length(goodLon) length(goodLat) 1 1])');
    end
end

%% subset and load precipitation forecasts

% only use 0Z and 12Z WPC forecasts
if initTime.Hour == 6 || initTime.Hour == 18
    initTimeWPC = initTime - hours(6);
else
    initTimeWPC = initTime;
end
timeP = (initTimeWPC + hours(incp)):hours(incp):(initTimeWPC + hours(lenp));

% get precipitation grid and subset indices
pfnametest = sprintf('%sp06m_%sf%03u.nc',precipPath,datestr(initTimeWPC,'yyyymmddhh'),incp);
% increment until finding a file that exists
i = incp;
while ~exist(pfnametest,'file')
    i = i + incp;
    pfnametest = sprintf('%sp06m_%sf%03u.nc',precipPath,datestr(initTimeWPC,'yyyymmddhh'),i);
    if i > ntP*incp
        disp('No vaild WPC precipitation files.');
        break
    end
end
lonPall = double(ncread(pfnametest,'gridlon_0'));
latPall = double(ncread(pfnametest,'gridlat_0'));
goodLonP = find(sum(lonPall >= boundingBox(1) & lonPall <= boundingBox(2),2) > 0);
goodLatP = find(sum(latPall >= boundingBox(3) & latPall <= boundingBox(4)) > 0);
lonWPC = lonPall(goodLonP,goodLatP);
latWPC = latPall(goodLonP,goodLatP);

% read in precip grids
precip = NaN(ntP,size(lonWPC,1),size(lonWPC,2));
for i = 1:ntP
    pfname = sprintf('%sp06m_%sf%03u.nc',precipPath,datestr(initTimeWPC,'yyyymmddhh'),incp*i);
    if exist(pfname,'file')
        if i == 1
            precip(i,:,:) = ncread(pfname,'APCP_P8_2L1_GLC0_acc',...
                [goodLonP(1) goodLatP(1)],[length(goodLonP) length(goodLatP)]);
        else
            precip(i,:,:) = ncread(pfname,'APCP_P8_2L1_GLC0_acc6h',...
                [goodLonP(1) goodLatP(1)],[length(goodLonP) length(goodLatP)]);
        end
    end
end

%% load and subset DEM

demObj = matfile(DEMfile);
latAll = demObj.Lat(:,1);
lonAll = demObj.Lon(1,:);
goodLon2 = find(lonAll > boundingBox(1) & lonAll < boundingBox(2));
goodLat2 = find(latAll > boundingBox(3) & latAll < boundingBox(4));
DEM = demObj.DEM(goodLat2,goodLon2);
DEM(isnan(DEM)) = 0;
[LonDEM,LatDEM] = meshgrid(lonAll(goodLon2),latAll(goodLat2));
clear demObj;

%% process West Coast HUC-8 watershed shapefiles into DEM masks

load(maskFile);
maskClipped = mask(:,:,basinsClippedInd);
basinsClipped = basins(basinsClippedInd);
basinsHypsoClipped = basinsHypso(:,basinsClippedInd);

% compute elevations at basin area deciles
basinDecileElm = NaN(length(basinsClipped),11);
basinsHypso2 = [zeros(1,length(basinsClipped)); basinsHypsoClipped];
bands2 = [0 bands];
ytick = 0.1:0.1:0.9;
for i = 1:length(basinsClipped)
    if ~isnan(basinsHypsoClipped(1,i))
        basinDecileElm(i,1) = bands(find(basinsHypsoClipped(:,i) > 0,1,'first')) - 50;
        bandsTemp = bands2;
        bandsTemp(find(roundn(basinsHypso2(:,i),-6) == 0,1,'last')) = bandsTemp(find(roundn(basinsHypso2(:,i),-6) == 0,1,'last')) + 50;
        basinDecileElm(i,2:10) = interp1(basinsHypso2(:,i) + cumsum((1e-6)*ones(46,1)),bandsTemp ,ytick);
        if basinsHypsoClipped(1,i) ~= 1
            basinDecileElm(i,11) = bands(find(roundn(basinsHypsoClipped(:,i),-6) < 1,1,'last')) + 50;
        else
            basinDecileElm(i,11) = bands(1);
        end
    end
end

%% bias-correct Z0C ensemble spread according to results in Henn et al. 2018

% establish bias-correction ratio as a function of time
%ratioCoeffs = [-0.8490 6.0875];
ratioCoeffs = [-0.4 4];
SDratio = polyval(ratioCoeffs, log1p(0:inct:lenf));

Z0C_min = 0;
Z0C_max = 5000;

Z0C_corrected = NaN(size(Z0C));
for i = 1:nt
    Z0C_mean = permute(repmat(squeeze(mean(Z0C(:,i,:,:))),1,1,nEns),[3 1 2]);
    Z0C_anom = squeeze(Z0C(:,i,:,:)) - Z0C_mean;
    Z0C_corrected(:,i,:,:) = max(min(Z0C_anom*SDratio(i) + Z0C_mean,Z0C_max),Z0C_min);
    fprintf('Finished with bias-correcting ensemble spread for %s.\n', time(i));
end

%% downscale freezing levels to fine DEM and compute basin values

probTerrainAboveZ0C = NaN(nt,size(LatDEM,1),size(LatDEM,2));
basinFracAboveZ0C = NaN(nt,nEns,length(basinsClipped));
basinMeanFracAboveZ0C = NaN(nt,length(basinsClipped));
basinZ0C = NaN(nt,nEns,length(basinsClipped));
basinMeanZ0C = NaN(nt,length(basinsClipped));
% pre-compute number of grid cells in each watershed for speed
nCellsBasin = NaN(length(basinsClipped),1);
for i = 1:length(basinsClipped)
    nCellsBasin(i) = sum(sum(maskClipped(:,:,i)));
end
maxCells = max(nCellsBasin);
DEMbasin = NaN(maxCells,length(basinsClipped));
for i = 1:length(basinsClipped)
   temp = DEM(maskClipped(:,:,i));
   DEMbasin(1:nCellsBasin(i),i) = temp;
end
for i = 1:nt
    Z0C_fine = NaN(size(LatDEM,1),size(LatDEM,2),nEns);
    terrainAboveZ0C = false(size(LatDEM,1),size(LatDEM,2),nEns);
    for j = 1:nEns
        temp = interp2(LonGEFS,LatGEFS,squeeze(Z0C_corrected(j,i,:,:)),LonDEM,LatDEM,'bilinear');
        Z0C_fine(:,:,j) = temp;
        terrainAboveZ0C(:,:,j) = temp - zOffset - DEM < 0;
        for k = 1:length(basinsClipped)
            temp3 = temp(maskClipped(:,:,k));
            basinZ0C(i,j,k) = mean(temp3);
            basinFracAboveZ0C(i,j,k) = sum(temp3 - zOffset - DEMbasin(1:nCellsBasin(k),k) < 0)/nCellsBasin(k);
        end
        fprintf('Finished with downscaling forecast ensemble member %u for timestep %u of %u.\n',j,i,nt);
    end
    for k = 1:length(basinsClipped)
        basinMeanZ0C(i,k) = mean(squeeze(basinZ0C(i,:,k)));
        basinMeanFracAboveZ0C(i,k) = mean(squeeze(basinFracAboveZ0C(i,:,k)));
    end
    probTerrainAboveZ0C(i,:,:) = mean(terrainAboveZ0C,3);
end

% sort basin ensembles for faster plotting later
basinFracAboveZ0Csorted = NaN(size(basinFracAboveZ0C));
for i = 1:length(basinsClipped)
    for j = 1:nt
        basinFracAboveZ0Csorted(j,:,i) = sort(basinFracAboveZ0C(j,:,i));
    end
end

%% downscale course precipitation data to basin averages at each 6-hrly timestep

basinPrecip = NaN(ntP,length(basinsClipped));
basinRain = NaN(ntP,length(basinsClipped));
basinSnow = NaN(ntP,length(basinsClipped));
basinMixd = NaN(ntP,length(basinsClipped));
precipFine = NaN(nt,size(LatDEM,1),size(LatDEM,2));
for i = 1:ntP
    F = scatteredInterpolant(lonWPC(:),latWPC(:),reshape(squeeze(precip(i,:,:)),length(lonWPC(:)),1),'linear');
    temp = F(LonDEM,LatDEM);
    precipFine(i,:,:) = temp;
    if initTime == initTimeWPC
        offset = 1;
    else
        offset = 0;
    end
    rain = zeros(size(temp));
    rain(squeeze(probTerrainAboveZ0C(i + offset,:,:)) == 0) = temp(squeeze(probTerrainAboveZ0C(i + offset,:,:)) == 0);
    snow = zeros(size(temp));
    snow(squeeze(probTerrainAboveZ0C(i + offset,:,:)) == 1) = temp(squeeze(probTerrainAboveZ0C(i + offset,:,:)) == 1);
    mixd = NaN(size(temp,1),size(temp,2),nEns - 1);
    for j = 1:(nEns - 1)
        temp2 = zeros(size(temp));
        temp2(squeeze(probTerrainAboveZ0C(i + offset,:,:)) >= j/nEns & squeeze(probTerrainAboveZ0C(i + offset,:,:)) < (j + 1)/nEns) = ...
            temp(squeeze(probTerrainAboveZ0C(i + offset,:,:)) >= j/nEns & squeeze(probTerrainAboveZ0C(i + offset,:,:)) < (j + 1)/nEns);
        mixd(:,:,j) = temp2;
    end
    for j = 1:length(basinsClipped)
        basinPrecip(i,j) = mean(temp(maskClipped(:,:,j)));
        basinRain(i,j) = mean(rain(maskClipped(:,:,j)));
        basinSnow(i,j) = mean(snow(maskClipped(:,:,j)));
        for k = 1:(nEns - 1)
            temp2 = squeeze(mixd(:,:,k));
            basinMixd(i,j,k) = mean(temp2(maskClipped(:,:,j)));
        end
        fprintf('Finished with downscaling precipitation forecast for basin %u of %u, timestep %u of %u.\n',...
            j,length(basinsClipped),i,ntP);
    end
end

%% create map forecast figures of ensemble probability of being below freezing at each time

% interpolate grids with scale factor to produce an "effective" Mercator
% projection for plotting on DS site

bottomInd = find((abs(LatDEM(:,1) - boundingBox(3))) == min(abs(LatDEM(:,1) - boundingBox(3))));
topInd = find((abs(LatDEM(:,1) - boundingBox(4))) == min(abs(LatDEM(:,1) - boundingBox(4))));
scaleFactor = cumsum(sec(deg2rad(LatDEM(bottomInd:-1:topInd,1))));
latMerc = (boundingBox(4) - boundingBox(3))*(scaleFactor - scaleFactor(1))/(scaleFactor(end) - scaleFactor(1)) + boundingBox(3);
LatDEMmerc = repmat(flipud(latMerc),1,size(LatDEM,2));
LonDEMmerc = LonDEM(topInd:bottomInd,:);

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultAxesFontSize',10);
load(cmapFile);
cmap2 = flipud(parula(10));
for i = 1:nt
    f = figure;
    f.Units = 'centimeters';
    f.Position = [2 2 22.1052 30];
    ax1 = axes;
    colormap(ax1,cmap);
    hold(ax1,'on');
    h = pcolor(LonDEMmerc,LatDEMmerc,squeeze(probTerrainAboveZ0C(i,topInd:bottomInd,:)));
    h.EdgeColor = 'none';
    ax2 = axes();
    colormap(ax2,cmap2);
    v = [1 5 12.5 25];
    if i > 1 && i <= (ntP + 1)
        [ch,hp] = contour(ax2,LonDEMmerc,LatDEMmerc,squeeze(precipFine(i - 1,topInd:bottomInd,:)),v);
        hp.LineWidth = 1;
        clabel(ch,hp,v,'Color',[0.6 0.6 0.6],'FontName','Arial');
    end
    ax2.CLim = [1 25];
    ax1.XTick = [];
    ax1.YTick = [];
    ax1.XLim = [boundingBox(1) boundingBox(2)];
    ax1.YLim = [boundingBox(3) boundingBox(4)];
    ax1.Position = [0 0 1 1];
    ax2.XLim = [boundingBox(1) boundingBox(2)];
    ax2.YLim = [boundingBox(3) boundingBox(4)];
    ax2.XTick = [];
    ax2.YTick = [];
    ax2.Position = ax1.Position;
    ax2.Color = 'none';
    %print(f,[archiveDir '/Map_' datestr(initTime,'yyyy-mm-dd_HH') '_'...
    %    datestr(initTime + hours(inct*(i - 1)),'yyyy-mm-dd_HH') '.png'],'-dpng','-r150');
    print(f,[outputDir '/Map_' sprintf('%03u',(inct*(i - 1))) 'hrs_current.png'],'-dpng','-r300');
    close(f);
end

%% create basin forecast figures of ensemble fraction of basin area below freezing over time 

dateticks = dateshift(time(1),'start','day'):days(1):dateshift(time(end),'start','day');
dateticklabels = cell(length(dateticks),1);
ylabels2 = {'0: Below Watershed','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1: Above Watershed'};
cmap3struct = load(cmap3file,'CoolWarmFloat257');
cmap3 = flipud(cmap3struct.CoolWarmFloat257(round(linspace(1,257,nEns)),:));
barColors = [0 0.5 0; 0.8 1 0.8; 1 1 1];
for i = 1:length(dateticks)
    dateticklabels{i} = sprintf('%s',datestr(dateticks(i),'dd-mmm'));
end
for i = 1:length(basinsClipped)
    basinDecileElLabels = cell(11,1);
    for j = 1:11
        basinDecileElLabels{j} = ...
            sprintf('%u m (%u ft)',round(basinDecileElm(i,j)),round(m2ft(basinDecileElm(i,j))));
    end 
    f2 = figure;
    f2.Units = 'centimeters';
    f2.Position = [2 2 25 15];
    ax1 = axes();
    hold(ax1,'on');
    plot(time,1 - squeeze(basinFracAboveZ0C(:,:,i)),'-','Color',[0.75 0.75 0.75]);
    % color-code ensemble members
    for j = 1:nEns
        he(j) = plot(time,1 - basinFracAboveZ0Csorted(:,j,i),'.');
        he(j).Color = cmap3(j,:);
        he(j).MarkerSize = 10;
    end
    hm = plot(time,1 - squeeze(basinMeanFracAboveZ0C(:,i)),'k-.','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','k');
    ax1.XTick = datenum(dateticks);
    ax1.XTickLabel = dateticklabels;
    ax1.XGrid = 'on';
    ax1.XLim = [datenum(initTime) datenum(initTime + hours(168))];
    ax1.YLim = [-0.01 1.01];
    ax1.YTick = 0:0.1:1;
    ax1.YTickLabel = basinDecileElLabels;
    ax1.YLabel.String = 'Forecast Freezing Level';
    ax1.XLabel.String = 'Date [UTC, PST + 8]';
    ax1.Title.String = sprintf('%s Forecast Initialized %sZ\n7-day WPC Precipitation Total: %4.1f mm (%4.2f in) - %u%% Rain, %u%% Rain or Snow, %u%% Snow',...
        basinsClipped(i).Name,datestr(initTime,'dd-mmm hh'),roundn(nansum(basinPrecip(:,i)),-1),...
        roundn(nansum(basinPrecip(:,i))/25.4,-2),round(100*nansum(basinRain(:,i))/nansum(basinPrecip(:,i))),...
        round(100*nansum(basinMixd(:,i))/nansum(basinPrecip(:,i))),round(100*nansum(basinSnow(:,i))/nansum(basinPrecip(:,i))));
    ax1.Title.FontSize = ax1.FontSize;
    ax2 = axes('Position',ax1.Position);
    ax2.Color = 'none';
    ax2.XTick = [];
    ax2.XLim = [datenum(initTime) datenum(initTime + hours(168))];
    ax2.YTick = 0:0.1:1;
    ax2.YTickLabel = ylabels2;
    ax2.YLabel.String = 'Fraction of Watershed Below Freezing Level';
    ax2.YLabel.Position(1) = datenum(time(29) + hours(13));
    ax2.YLabel.Rotation = 270;
    ax2.YLim = [-0.01 1.01];
    ax2.YAxisLocation = 'right';
    ax1.Position = [0.20 0.15 0.52 0.75];
    ax2.Position = ax1.Position;
    ax3 = axes;
    hold(ax3,'on');
    y3max = 1.3*nansum(basinPrecip(:,i));
    hb = bar(datenum(timeP - hours(incp/2)),[basinRain(:,i) squeeze(basinMixd(:,i,:)) basinSnow(:,i)],'stacked');
    for j = 1:(nEns + 1)
        hb(j).FaceColor = barColors(j,:);
        hb(j).EdgeColor = 'none';
    end
    hb2 = bar(datenum(timeP - hours(incp/2)),basinPrecip(:,i));
    hb2.FaceColor = 'none';
    hb2.EdgeColor = 'k';
    for j = 1:ntP
        if roundn(basinPrecip(j,i),-1) > 0.254
            ht = text(datenum(timeP(j)) - 7.5/24,basinPrecip(j,i) + 0.03*y3max,...
                sprintf('%4.1f',roundn(basinPrecip(j,i),-1)));
            ht.FontSize = 10;
            ht.Color = [0 0 0];
        end
    end
%     htt = text(datenum(timeP(end)) - 4.5,0.95*y3max,...
%         sprintf('8-day GEFS median precipitation total: %4.1f mm (%4.2f in)',...
%         roundn(sum(basinMedianPrecip(:,i)),-1),roundn(sum(basinMedianPrecip(:,i))/25.4,-2)));
%     htt.FontSize = 10;
%     htt.Color = [0 0 1];
%     hpa = plot(timeP + hours(3),cumsum(basinPrecip(:,i)),'--');
%     hpa.LineWidth = 2;
%     hpa.Color = [0 0 1];
    ax3.XLim = [datenum(initTime) datenum(initTime + hours(168))];
    ax3.YLim = [0 y3max + 1e-6];
    ax3.XTick = [];
    ax3.YTick = [];
    ax3.Color = 'none';
    ax3.Position = ax1.Position;
    hl = legend(ax3, [hm he(1) he(nEns) hb(1) hb(ceil(nEns/2)) hb(nEns + 1)],'GEFS Mean','GEFS Warmest','GEFS Coolest',...
        'All Rain','Half Rain, Half Snow','All Snow');
    hl.FontSize = ax1.FontSize;
    hl.Orientation = 'horizontal';
    hl.Position = [0.1829 0.0347 0.5450 0.0344];
    % add second plot showing hypsometric curve
    ax4 = axes;
    hold(ax4,'on');
    hh = plot(0:0.1:1,basinDecileElm(i,:),'k-');
    ax4.XGrid = 'on';
    ax4.XLim = [0 1];
    ax4.YLim = [0 4100];
    ax4.XTick = 0:0.25:1;
    ax4.YTick = 0:500:4000;
    ax4.XLabel.String = 'Fraction of Watershed Below';
    ax4.YLabel.String = 'El. [m]';
    ax4.Title.String = sprintf('%s Hypsometry',basinsClipped(i).Name);
    ax4.Position = [0.82 0.25 0.12 0.55];
    ax4.Box = 'off';
    ax5 = axes;
    ax5.Color = 'none';
    ax5.XTick = [];
    ax5.XLim = [0 1];
    ax5.YTick = 0:2000:12000;
    ax5.YLabel.String = 'El. [ft]';
    ax5.YLabel.Rotation = 270;
    ax5.YLabel.Position(1) = 1.4;
    ax5.YLabel.Position(2) = 6500;
    ax5.YLim = 3.28*[0 4100];
    ax5.YAxisLocation = 'right';
    ax5.Box = 'off';
    ax5.Position = ax4.Position;
    ax6 = axes;
    ax6.Color = 'none';
    ax6.XTick = [];
    ax6.YTick = [];
    ax6.Position = ax4.Position;
    %print(f2,[archiveDir '/' num2str(basinsClipped(i).HUC8) '_' datestr(initTime,'yyyy-mm-dd_HH') '.png'],'-dpng','-r150');
    print(f2,[outputDir '/' num2str(basinsClipped(i).HUC8) '_current.png'],'-dpng','-r300');
    close(f2);
end
clear he;

%%

toc;
quit;
