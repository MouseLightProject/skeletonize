function varargout = clusterSkelGenPipeline(inputFile, pipelineInputFile)
% CLUSTERSKELGENPIPELINE Create Mouse Light Pipeline input file.
%
%   CLUSTERSKELGENPIPELINE(INPUTFILE, DATASET, PIPELINEINPUTPUTFILE) creates a Mouse Light Pipeline input file for the
%   given HDF5 source file.
%
%   This is an alternative to the portion of the original cluster_skelh5.m that generated a batch file for submitting
%   all jobs to the cluster.  It generates the required Mouse Light Pipeline input file for processing.
%
%   INPUT PARAMETERS:
%       INPUTFILE:            absolute path to configuration file (see notes).
%       DATASET:              dataset name in source file.
%       PIPELINEINPUTPUTFILE: output Pipeline input file name.
%
%   It assumed the support files in the skeletonization repository (./common and downstream) are on the MATLAB path.
%
%   Notes:
%       The path to the configuration file and the path specified in the configuration file to the HDF5 data file will
%       be passed from the pipeline to the skeletonization task.  If the path on the system used to generate this file
%       is different than the pipeline system has access to (e.g., this script is used on Windows and accesses nrs and
%       dm11 via a mapped drive or \\dm11... and \\nrs vs. /nrs/... and /groups/... as seen by the pipeline) you will
%       need to modify the paths in the customParameters entry in the generated JSON file.
%
%   Examples:
%       clusterSkelGenPipeline('/groups/mouselight/mouselight/pipeline-systems/support/skeletonization/config_files/20170925_prob0_config_skelh5.cfg', 'pipeline-input.json');

opt = configparser(inputFile);

[brainSize,RR,chunk_dims] = h5parser(opt.inputh5, opt.h5prob);

% get a multiple of chunksize that is around 1000^3
cropSize = round(1000./chunk_dims).*chunk_dims;

% to get %10 overlap overhead use multiple of 10
fullh = chunk_dims; % add 1 to make it odd (heuristic)

outfolder = fileparts(pipelineInputFile);

if numel(outfolder) > 0 &&  ~exist(outfolder, 'dir')
    mkdir(outfolder);
end

bbox = createOverlapBox(brainSize,cropSize,fullh);

BBoxes = RR(:,[1 4 2 5 3 6])+1;

X = BBoxes(:,1:2);
Y = BBoxes(:,3:4);
Z = BBoxes(:,5:6);
XYZ = unique([X(:),Y(:),Z(:)],'rows');
in = inhull([bbox(:,1:2:end);bbox(:,2:2:end)],XYZ);
in = any(reshape(in,[],2),2);

numTilesMax = size(bbox,1);

fid = fopen(pipelineInputFile,'w');

projectFormat = struct('major', 1, 'minor', 0);

ts = struct('id', num2cell(1:numTilesMax), ...
    'relativePath', arrayfun(@(x) sprintf('%05d', x), 1:numTilesMax, 'UniformOutput', false), ...
    'isComplete', true, ...
    'position', struct('x', {}, 'y', {}, 'z', {}),...
    'step', struct('x', {}, 'y', {}, 'z', {}));

for idx = 1:numTilesMax
    BB = bbox(idx,:);
    
    % check if BB is outsize of BBoxes
    if ~in(idx)
        continue
    end
    
    ts(idx).position = struct('x', BB(1), 'y', BB(3), 'z', BB(5));
    ts(idx).step = struct('x', BB(2), 'y', BB(4), 'z', BB(6));
end

ts = ts(~cellfun(@isempty, {ts.position}));

position = [ts.position];
x = [position.x];
y = [position.y];
z = [position.z];

projectInfo = struct('customParameters', struct('inputFile', strrep(opt.inputh5, '\', '\\'), 'dataset', opt.h5prob, 'configurationFile',  strrep(inputFile, '\', '\\')));

content = struct('pipelineFormat', projectFormat, 'projectInfo', projectInfo,  ...
    'extents', struct('minimumX', min(x), 'maximumX', max(x), 'minimumY', min(y), 'maximumY', max(y),'minimumZ', min(z), 'maximumZ', max(z)), ...
    'tiles', ts);

fprintf(fid, jsonencode(content));

fclose(fid);

if nargout > 0
    varargout{1} = content;
end
