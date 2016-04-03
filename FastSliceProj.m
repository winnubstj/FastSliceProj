function [ ] = FastSliceProj(  )
%FastSliceProj Takes JSON file from Acquisition dashboard and generates projection off all stacks belonging to same slice.
% Dependency: fastProtoBuf

tileDir = 'Y:\mousebrainmicro\acquisition\3-21-2016\Data\';
reqZ = 172;

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
fprintf('Done\n');

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
FOV.x_res = FOV.x_size_um/FOV.x_size_pix;
FOV.y_res = FOV.y_size_um/FOV.y_size_pix;


%% Get general dimensions result file.
fileLoc(latPos(:,3)==reqZ,:);
latPos(latPos(:,3)==reqZ,:);

end

