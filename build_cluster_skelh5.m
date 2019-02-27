if ~exist(fullfile(pwd(), 'compiled'), 'file') ,
    mkdir(fullfile(pwd(), 'compiled')) ;
end
mcc -m -R -nojvm -v cluster_skelh5.m -d ./compiled  -a ./common
copyfile run_cluster_skelh5_at_janelia.sh ./compiled
