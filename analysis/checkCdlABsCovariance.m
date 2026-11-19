function summary = checkCdlABsCovariance(outputDir)
%CHECKCDLABSCOVARIANCE Check Monte Carlo variation in CDL-A receive correlation.
% The 9x9 panel retains every port pair used by the 16-port covariance metric.
arguments
    outputDir (1,1) string = fullfile(fileparts(mfilename('fullpath')),"cdl_toolbox_comparison")
end
savedRng = rng;
restoreRng = onCleanup(@()rng(savedRng));
cfg = configUL4T('base');
cfg.Nh = 9;
cfg.Nv = 9;
cfg.users.mode = 'specified';
cfg.users.positions_m = [0 60 1.5];
cfg.ul.ue_speed_kmh = 30;
cfg.cdl.toolbox_profile = 'CDL-A';
cfg = prepareFixedArrayConfig(cfg);
pattern = fixedArrayChannel('pattern',cfg);
schedule = struct('srs_ms',0,'ul_ms',(0:1:40).','slot_ms',1);
ue = struct('layout','upa','positions_lambda',.22*[0 0 0;0 0 1;1 0 0;1 0 1], ...
    'antenna_indices',1:4,'toolbox_ports',1:4);
[row,col] = ndgrid([1 2 5 9],[1 2 5 9]);
ports = sub2ind([cfg.Nv cfg.Nh],row(:),col(:));

bankSeeds = [26101001 26201001];
subsetSizes = [40;80;160];
covSum = zeros(numel(ports),numel(ports),numel(bankSeeds),2,numel(subsetSizes));
for bank = 1:numel(bankSeeds)
    cfg.seed = bankSeeds(bank);
    accumulating = zeros(numel(ports),numel(ports),2);
    for mc = 1:max(subsetSizes)
        rng(cfg.seed+mc-1,'twister');
        native = ul4tNativeCDL(cfg,pattern,schedule,mc,ue);
        rng(cfg.seed+mc-1,'twister');
        ref = ul4tToolboxChannel(cfg,pattern,schedule,mc,ue);
        accumulating(:,:,1) = accumulating(:,:,1)+receiveGram(native,ports);
        accumulating(:,:,2) = accumulating(:,:,2)+receiveGram(ref,ports);
        cutoff = find(subsetSizes==mc,1);
        if ~isempty(cutoff)
            covSum(:,:,bank,:,cutoff) = reshape(accumulating,numel(ports),numel(ports),1,2);
        end
    end
    fprintf('CDL-A BS covariance bank %d finished, %d seeds\n',bank,max(subsetSizes));
end

summary = table();
for i = 1:numel(subsetSizes)
    a = covSum(:,:,:,1,i);
    b = covSum(:,:,:,2,i);
    summary = [summary;table(subsetSizes(i),2*subsetSizes(i), ...
        100*distance(a(:,:,1),b(:,:,1)),100*distance(a(:,:,2),b(:,:,2)), ...
        100*distance(a(:,:,1),a(:,:,2)),100*distance(b(:,:,1),b(:,:,2)), ...
        100*distance(sum(a,3),sum(b,3)), ...
        'VariableNames',{'SeedsPerBank','SeedsPooled','BaseVsToolboxBank1_pct', ...
        'BaseVsToolboxBank2_pct','BaseAcrossBanks_pct','ToolboxAcrossBanks_pct', ...
        'BaseVsToolboxPooled_pct'})]; %#ok<AGROW>
end
if ~isfolder(outputDir), mkdir(outputDir); end
writetable(summary,fullfile(outputDir,'cdl_a_bs_covariance_convergence.csv'));
save(fullfile(outputDir,'cdl_a_bs_covariance_convergence.mat'),'summary','covSum','bankSeeds');
disp(summary);
end

function G = receiveGram(model,ports)
gains = model.pathGains(:,:,:,ports);
pattern = single(sqrt(model.diagnostics.cluster_pattern_power));
data = reshape(gains.*reshape(pattern,1,[]),[],numel(ports));
G = double(data'*data);
end

function d = distance(G,H)
% Symmetric relative difference so that exchanging the banks has no effect.
G = G./sqrt(real(diag(G))*real(diag(G)).');
H = H./sqrt(real(diag(H))*real(diag(H)).');
d = 2*norm(G-H,'fro')/(norm(G,'fro')+norm(H,'fro'));
end
