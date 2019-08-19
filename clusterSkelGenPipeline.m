function clusterSkelGenPipeline(inputFile, dataset, pipelineInputFile)
% CLUSTERSKELGENPIPELINE Create Mouse Light Pipeline input file.
%
%   CLUSTERSKELGENPIPELINE(INPUTFILE, DATASET, PIPELINEINPUTPUTFILE) creates
%   a Mouse Light Pipeline input file for the given HDF5 source file.
%
%   This is an alternative to the portion of the original cluster_skelh5.m
%   that generated a batch file for submitting all jobs to the cluster.  It
%   generates the required Mouse Light Pipeline input file for processing.
%
%   INPUT PARAMETERS:
%       INPUTFILE:            source HDF5 file.
%       DATASET:              dataset name in source file.
%       PIPELINEINPUTPUTFILE: output Pipeline input file name.
%
%   Examples:
%       clusterSkelGenPipeline('/nrs/mouselight/cluster/classifierOutputs/2017-09-25/20170925_prob0/20170925_prob0_lev-6_chunk-111_111_masked-0.h5',
%       '/prob0', 'pipeline-input.json');

[brainSize,RR,chunk_dims] = h5parser(inputFile,dataset);

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

pf = struct('major', 1, 'minor', 0);
pi = struct('inputFile', strrep(inputFile, '\', '\\'), 'dataset', dataset);

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

fprintf(fid, jsonencode(struct('pipelineFormat', pf, 'projectInfo', pi, 'tiles', ts)));

fclose(fid);
