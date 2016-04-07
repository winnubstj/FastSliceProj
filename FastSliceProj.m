function [ ] = FastSliceProj(  )
%FastSliceProj Takes JSON file from Acquisition dashboard and generates projection off all stacks belonging to same slice.
% Dependency: fastProtoBuf, readTifFast, writeTifFast

tileDir = 'Y:\mousebrainmicro\acquisition\3-21-2016\Data\';
tempDir = 'C:\Users\winnubstj\Desktop\Temp';
channels = [0,1];
reqZ = 187;
frames = 100:150;
downSample = 2;

%% load and convert Json file to .mat
fprintf('\n');
jsonFile = fullfile(tileDir,'dashboard.json');
matFile = fullfile(tileDir,'dashboard.mat');
if isempty(dir(matFile))
    %load JSON
    fprintf('Loading Json file..\n');
    if isempty(dir(jsonFile)), error('JSON file: %s could not be found',jsonFile); end
    opt.ShowProgress = true;
    jsonData = loadjson(jsonFile,opt);
    fprintf('Saving Json file..\n');
    save(matFile,'jsonData');
else
    %load Mat
    fprintf('Loading pre-parsed Json file..\n');
    load(matFile);
end
fprintf('Done!\n');

%% Get lattice position and file location all tiles.
contents = structfun(@(x) x,jsonData.tileMap,'UniformOutput',false);
contents = struct2cell(contents);
nTiles = cellfun(@(x) size(x,2),contents); %only days that have tiles.
contents = contents(nTiles>1);
contents = horzcat(contents{:});
latPos = cellfun(@(x) x.contents.latticePosition,contents); % get Lattice position.
latPos = [latPos.x;latPos.y;latPos.z]';
fileLoc = cellfun(@(x) fullfile(jsonData.monitor.location,x.relativePath),contents,'UniformOutput',false)'; % get file position.

%% Get rid of duplicates (use most recent)
latMat = zeros(max(latPos(:,1)),max(latPos(:,2)),max(latPos(:,3)));
for iTile = 1:size(latPos,1)
    if latPos(iTile,1)>0
        latMat(latPos(iTile,1),latPos(iTile,2),latPos(iTile,3)) = latMat(latPos(iTile,1),latPos(iTile,2),latPos(iTile,3)) +1;
    end
end
% find lattice positions with duplicate entries.
[x,y,z] = ind2sub(size(latMat),find(latMat>1));removeList = [];
% Select earliest tile for removal.
for iTile = 1:size(x,1)
   ind = find(latPos(:,1)==x(iTile) &  latPos(:,2)==y(iTile) &  latPos(:,3)==z(iTile));
   removeList = [removeList; min(ind)];
end
fileLoc(removeList,:) = [];
latPos(removeList,:) = [];

%% Get field of view info from protobuf file.
ind = find(latPos(:,3)>0,1);
FOV = fastProtoBuf( fullfile(fileLoc{ind},[fileLoc{ind}(end-4:end),'-ngc.microscope']), {'x_size_um','y_size_um','z_size_um','x_overlap_um','y_overlap_um','z_overlap_um'});
FOV.x_size_pix = contents{ind}.contents.imageContents{1,1}.width;
FOV.y_size_pix = contents{ind}.contents.imageContents{1,1}.height;
FOV.z_size_pix = contents{ind}.contents.imageContents{1,1}.frameCount;
FOV.x_res = FOV.x_size_um/FOV.x_size_pix;
FOV.y_res = FOV.y_size_um/FOV.y_size_pix;
FOV.x_overlap_pix = FOV.x_overlap_um/FOV.x_res;
FOV.y_overlap_pix = FOV.y_overlap_um/FOV.y_res;

%% Get general dimensions result file.
fileLoc = fileLoc(latPos(:,3)==reqZ,:);
latPos = latPos(latPos(:,3)==reqZ,:);
latPos(:,1) = latPos(:,1)-min(latPos(:,1))+1;
latPos(:,2) = latPos(:,2)-min(latPos(:,2))+1;
tileWidth = max(latPos(:,2)) - min(latPos(:,2));
tileHeight = max(latPos(:,1)) - min(latPos(:,1));
fprintf('Proj. will consist of %i tiles (%i width and %i height)\n',length(fileLoc),tileWidth,tileHeight);
% Display tile map so user can check for missing.
latMat = zeros(max(latPos(:,1)),max(latPos(:,2)));
latMat(sub2ind(size(latMat),latPos(:,1),latPos(:,2)))=1;
figure(); imshow(latMat,[]);

%% Create Projections
if isempty(dir(tempDir)),mkdir(tempDir); end
fprintf('Starting Projection cycle\n');
nTiles = size(fileLoc,1);
% nTiles = 73;
IStack = zeros([FOV.y_size_pix,FOV.x_size_pix,nTiles,length(channels)],'uint16');
mmPos = [];
for iTile = 1:nTiles
    cTileDir = fullfile(tempDir,fileLoc{iTile}(end-18:end));
    if isempty(dir(cTileDir)),mkdir(cTileDir); end
    pos = fastProtoBuf( fullfile(fileLoc{iTile},[fileLoc{iTile}(end-4:end),'-ngc.acquisition']), {'x_mm','y_mm','z_mm'});
    mmPos = [mmPos;pos.x_mm,pos.y_mm,pos.z_mm];
    % Check both channels.
    for iChan = channels
       fprintf('%s - [%i/%i Chan: %i] ',datestr(now,'HH:MM'), iTile,size(fileLoc,1),iChan);
       %check if file already excists.
       cTileSource = fullfile(fileLoc{iTile},[fileLoc{iTile}(end-4:end),'-ngc.',num2str(iChan),'.tif']);
       cTileDest = fullfile(cTileDir,[cTileDir(end-4:end),'-ngc.',num2str(iChan),'.tif']);
       if isempty(dir(cTileDest))
           fprintf('Creating Projection\n');
           I = readTifFast(cTileSource,[FOV.y_size_pix,FOV.x_size_pix],frames,'uint16');
           I = max(I,[],3);
           writeTifFast(I,cTileDest,'uint16');
       else
           fprintf('Reading Stored Projection\n');
           I = readTifFast(cTileDest,[FOV.y_size_pix,FOV.x_size_pix],1,'uint16');
       end
       IStack(:,:,iTile,iChan+1) = I;
    end
end

%% Uneven illumination correction.
BG = squeeze(nanmedian(IStack(:,:,:,1),3));
BG = double(BG);
BG = BG./repmat(max(max(BG)),size(BG,1),size(BG,2));
BG = repmat(1,size(BG,1),size(BG,2))./BG;
% figure(); imshow(BG,[]);
% Temp nutil I put in overlap correct.
BG(BG>1.15)=1;
%% Blend tiles into full image
% Create spatial referencing object.
XExtentmm = [min(mmPos(:,1)),max(mmPos(:,1))+(FOV.x_size_um*10^-3) ];
YExtentmm = [min(mmPos(:,2)),max(mmPos(:,2))+(FOV.y_size_um*10^-3)];
XExtentPix = round((((XExtentmm(2)-XExtentmm(1))*1000)/FOV.x_res)/downSample);
YExtentPix = round((((YExtentmm(2)-YExtentmm(1))*1000)/FOV.y_res)/downSample);
resI = zeros(XExtentPix,YExtentPix,'uint16');
ROUT = imref2d(size(resI));
ROUT.XWorldLimits = XExtentmm;
ROUT.YWorldLimits = YExtentmm;

for iChan = 1:length(channels)
    resI = zeros(XExtentPix,YExtentPix,'uint16');
    for iTile = 1:nTiles
        tform = affine2d([1 0 0; 0 1 0; 0 0 1]);% Can add stage displacement object here
        RITile = imref2d([size(IStack,1), size(IStack,2)],[mmPos(iTile,1), mmPos(iTile,1)+(FOV.x_size_um*10^-3)],[mmPos(iTile,2), mmPos(iTile,2)+(FOV.y_size_um*10^-3)]) ; 
        IC = fliplr(flipud(IStack(:,:,iTile,iChan)));
        IC =double (IC);
        % Block overlap. 
        IC = uint16(IC.*BG);
        [J,RJ] = imwarp(IC,RITile,tform,'OutputView',ROUT);
        resI = max(cat(3,resI,J),[],3);
        iTile
    end
    imwrite(resI,fullfile(tempDir,['res.',num2str(iChan),'.tif']));
end

figure(); imshow(resI,[]);
imwrite(resI,fullfile(tempDir,'res.tif'));
end

