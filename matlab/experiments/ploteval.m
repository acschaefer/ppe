function ploteval
% PLOTEVAL Visualize results of plane extraction experiments.

%% Load data.
% Load evaluation data.
evaluation = load(fullfile('result','eval.mat'));

ngtall = evaluation.ngtall;
nmsall = evaluation.nmsall;
ncorrseg = evaluation.ncorrseg;
noverseg = evaluation.noverseg;
nundrseg = evaluation.nundrseg;
nmissseg = evaluation.nmissseg;
nnoisseg = evaluation.nnoisseg;
angdiff = evaluation.angdiff;
nangdiff = evaluation.nangdiff;
fcorrseg = evaluation.fcorrseg;
rmsems = evaluation.rmsems;
rmsegt = evaluation.rmsegt;
algorithmname = evaluation.algorithmname;

try
    e = evaluation.e;
    x = sqrt(e);
    xtitle = '$\sqrt{e}$';
end
try
    tcomps = evaluation.tcomps;
    x = tcomps;
    xtitle = '$T$';
end

if false
    rratio = rmsems(1,:,1,1) ./ rmsegt(1,:,1)
    ndiff = nmsall(1,:,1,1)-ngtall
    rratio = rratio(ndiff==0)
    rmean = mean(rratio)
    rstd = std(rratio)/sqrt(numel(rratio))
    [~,rt] = ttest(rratio-1)
end

%% Plot data.
makeplt('Number of reconstructed planes', 'Planes [1]', xtitle, ...
    'nms_real', x, ngtall, algorithmname, nmsall, 0);
makeplt('Correctly segmented on real data', 'Correctly segmented [1]', ...
    xtitle, 'c_real', x, ngtall, algorithmname, ncorrseg, 0);
makeplt('Correctly segmented on real data relative', ...
    'Correctly segmented $[1]$', xtitle, 'c_real_rel', x, ngtall, ...
    algorithmname, ncorrseg, 1);
makeplt('Over-segmented on real data', 'Over-segmented [1]', xtitle, ...
    'o_real', x, ngtall, algorithmname, noverseg, 0);
makeplt('Under-segmented on real data', 'Under-segmented [1]', xtitle, ...
    'u_real', x, ngtall, algorithmname, nundrseg, 0);
makeplt('Missing on real data', 'Missing $[1]$', xtitle, ...
    'm_real', x, ngtall, algorithmname, nmissseg, 0);
makeplt('Noise on real data', 'Noise $[1]$', xtitle, ...
    'n_real', x, ngtall, algorithmname, nnoisseg, 0);
makeplt('Angular difference', 'Angle [degree]', xtitle, ...
    'a_real', x, nangdiff, algorithmname, angdiff, 2);
makeplt('Fraction correctly segmented points', 'Correctly segmented [1]', ...
    xtitle, 'f_real', x, ngtall, algorithmname, fcorrseg, 0);
makeplt('RMSE reconstructed', 'RMSE $[\textup{m}]$', xtitle, ...
    'ems_real', x, ngtall, algorithmname, rmsems, 0);

if false
    ids = 1;
    rmsegt_real = sum(rmsegt(ids,:,1),2)/size(ngtall, 2);
    line([0,1],[rmsegt_real,rmsegt_real])
    fprintf('RMSE ground-truth: %f\n', rmsegt_real);
end

end

function makeplot(title, ytitle, xtitle, fname, x, ngtall, ...
    algorithmname, data)
    makeplt([title, ' relative'], ytitle, xtitle, [fname, '_rel'], x, ...
        ngtall, algorithmname, data, 1);
    makeplt([title, ' total'], ytitle, xtitle, [fname, '_tot'], x, ...
        ngtall, algorithmname, data, 0);
end


function makeplt(title, ytitle, xtitle, fname, x, denom, ...
    algorithmname, data, idenom)
    % Define line width in plots.
    lw = 1;

    % Plot data for each method
    fig = figure('Name', title);
    hold on
    
    ia = 1; % plane extraction method

    for ids = 2:size(data,1)
        % Denominator = number of point-clouds
        d = size(denom, 2);
        if idenom == 1
            % Denominator = total number of ground-truth planes
            d = sum(denom(ids,:));
        elseif idenom == 2
            % Denominator = total number of plane angle differences
            d = sum(denom(ids,:,ia,:),2);
        end

        dclean = data(ids,:,ia,:);
        if strcmp(fname, 'ems_real')
            dirty = dclean>0.02;
            ndirty = sum(dirty(:,:,:,:), 2);
            if max(ndirty) > 0
                warning('Removing <= %i outlier from %d total.', ...
                    max(ndirty), d);
                d = d-ndirty;
                dclean(dirty) = 0;
            end
        end
        y = squeeze(sum(dclean,2)./d);
        y(y==0) = NaN;

        plot(x, y, 'LineWidth',lw)

        ebest = 0.32; % SegComp
        if ids > 1
            ebest = 0.24; % Synthetic
        end
        ibest = find(abs(x-ebest)<1e-10);
        if numel(ibest) == 1
            fprintf('%s %i %f %f\n', fname, ids, x(ibest), y(ibest));
        end

    end
    hold off
    set(gca, 'TickLabelInterpreter', 'latex')
    xlabel(xtitle, 'Interpreter', 'latex')
    ylabel(ytitle, 'Interpreter', 'latex')
    legend(["$\phantom{0}0\ \textup{mm}$", ...
        "$\phantom{0}5\ \textup{mm}$", "$10\ \textup{mm}$", ...
        "$20\ \textup{mm}$", "$40\ \textup{mm}$"], ...
        'Interpreter', 'latex', 'Location', 'southeast')
    fname = fullfile('result',fname);
    savefig(fig, [fname, '.fig'])
    saveas(fig, [fname, '.png']);
    %matlab2tikz(sprintf('%s.tex', fname))
    
end
