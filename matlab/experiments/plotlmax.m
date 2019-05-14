function ploteval
% PLOTEVAL Visualize results of line extraction lmax scan.

%lmaxs = 0.006:0.002:0.02; % SegComp
lmaxs = 0.04:0.02:0.2; % Sim

cmax = [];
cbest = [];
obest = [];
ubest = [];
mbest = [];
nbest = [];
abest = [];
rbest = [];

for ilmax = 1:numel(lmaxs)
    lmax = lmaxs(ilmax);
    fprintf("lmax = %f\n", lmax);

    %% Load data.
    % Load evaluation data.
    %evaluation = load(fullfile('pcd','result','planeextract','train', ...
    %    'v1',['eval_lmax', num2str(lmax), '.mat']));
    evaluation = load(fullfile('pcd','result','planeextract','train_syn', ...
        ['eval_lmax', num2str(lmax), '.mat']));

    ngtall = evaluation.ngtall;
    nmsall = evaluation.nmsall;
    ncorrseg = evaluation.ncorrseg;
    noverseg = evaluation.noverseg;
    nundrseg = evaluation.nundrseg;
    nmissseg = evaluation.nmissseg;
    nnoisseg = evaluation.nnoisseg;
    angdiff = evaluation.angdiff;
    nangdiff = evaluation.nangdiff;
    rmsems = evaluation.rmsems;
    rmsegt = evaluation.rmsegt;
    algorithmname = evaluation.algorithmname;

    e = evaluation.e;
    x = sqrt(e);

    ids = 5;
    ngtall = squeeze(ngtall(ids,:));

    %% Plot data.
    [cmax, cbest] = getybest(x, ngtall, ncorrseg(ids,:,:,:), 1, cmax, cbest);
    [~, obest] = getybest(x, ngtall, noverseg(ids,:,:,:), 0, [], obest);
    [~, ubest] = getybest(x, ngtall, nundrseg(ids,:,:,:), 0, [], ubest);
    [~, mbest] = getybest(x, ngtall, nmissseg(ids,:,:,:), 0, [], mbest);
    [~, nbest] = getybest(x, ngtall, nnoisseg(ids,:,:,:), 0, [], nbest);
    [~, abest] = getybest(x, nangdiff(ids,:,:,:),  angdiff(ids,:,:,:), 2, [], abest);
    [~, rbest] = getybest(x, ngtall, rmsems(ids,:,:,:), 0, [], rbest);

end

lw = 1;
figure; plot(lmaxs, cmax, 'LineWidth', lw); ylabel('cmax');
figure; plot(lmaxs, cbest, 'LineWidth', lw); ylabel('cbest');
figure; plot(lmaxs, obest, 'LineWidth', lw); ylabel('obest');
figure; plot(lmaxs, ubest, 'LineWidth', lw); ylabel('ubest');
figure; plot(lmaxs, mbest, 'LineWidth', lw); ylabel('mbest');
figure; plot(lmaxs, nbest, 'LineWidth', lw); ylabel('nbest');
figure; plot(lmaxs, abest, 'LineWidth', lw); ylabel('abest');
figure; plot(lmaxs, rbest, 'LineWidth', lw); ylabel('rbest');

end


function [ymax, ybest] = getybest(x, ngtall, data, idenom, ymax, ybest)

    % Denominator = number of point-clouds
    denom = size(ngtall, 2);
    if idenom == 1
        % Denominator = total number of ground-truth planes
        denom = sum(ngtall);
    elseif idenom == 2
        % Denominator = total number of plane angle differences
        denom = sum(ngtall,2);
    end
    y = sum(data,2)./denom;

    % Best value is determined by specific e value
    %ebest = 0.32; % SegComp
    ebest = 0.24; % Sim
    ibest = find(abs(x-ebest)<1e-10);
    
    % Get maximum and best values
    ia = 1;
    ymax(end+1) = max(y(1,:,ia,:));
    ybest(end+1) = y(1,:,ia,ibest);
end
