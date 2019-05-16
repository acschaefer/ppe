function evalexp
% EVALEXP Evaluate plane extraction experiments.

%% Evaluate experimental results.
% Load dataset.
dataset = load(fullfile('..', 'data','dataset_test.mat'));
pc = dataset.pc;
gtang = dataset.gtang;

% Load first experimental result file.
plndata = load(fullfile('result','extrpln_ids1.mat'));
algorithmname = plndata.algorithmname;
pln = plndata.pln;
spln = size(pln);

% Initialize progress display.
it = 0;
nit = prod(spln(1:3));
disp('Evaluating results ...')

% Compare tolerance
tcomp = 0.8;

% Results are determined as function of the error threshold
e = (0.02:0.02:1).^2;

% Or fixed to #reconstructed = #ground-truth planes
fixngt = false;

% Scan compare tolerance (and fix e)
scantcomp = false;

% Plot the fits used to determine ground-truth RMSE
plotfits = false;

% Speed-up for fixngt
if fixngt
    e = [0];
end

% Initialize tcomp scan
if scantcomp
    tcomps = 0.51:0.01:1;
    ebest = 0.32^2;
    assert(isequal(size(tcomps), size(e)), ...
        'Array tcomp and e are not of the same size!');
    assert(~fixngt, 'fixngt and scantcomp both true!');
end

% Initialize results
ngtall = zeros(spln(1:2));
nmsall = zeros([spln(1:3),size(e,2)]);
ncorrseg = zeros([spln(1:3),size(e,2)]);
noverseg = zeros([spln(1:3),size(e,2)]);
nundrseg = zeros([spln(1:3),size(e,2)]);
nmissseg = zeros([spln(1:3),size(e,2)]);
nnoisseg = zeros([spln(1:3),size(e,2)]);
angdiff = zeros([spln(1:3),size(e,2)]);
nangdiff = zeros([spln(1:3),size(e,2)]);
fcorrseg = zeros([spln(1:3),size(e,2)]);
rmsems = zeros([spln(1:3),size(e,2)]);
rmsegt = zeros(spln(1:2));

% Loop over all datasets.
idss = 1:6;
for ids = idss
    fprintf("Running ids = %i\n", ids);

    % Load experimental results.
    plndata = load(fullfile('result',...
        ['extrpln_ids', num2str(ids), '.mat']));
    pln = plndata.pln;

    % Loop over all pointclouds.
    for ipc = 1 : size(pln,2)
        % Determine ground-truth indices
        pci = pc{ids, ipc};
        gti = unique(pci.Intensity);
        gti = gti(gti >= 10);
        ngt = size(gti, 1);
        ngtall(ids, ipc) = ngt;
        locall = reshape(pci.Location, [numel(pci.Location)/3, 3]);
        % Loop over all plane extraction methods.
        for ia = 1 : size(pln,3)
            % Loop over all parameters.
            ipmprev = -1;
            for ie = 1 : numel(e)
                ei = e(ie);
                if scantcomp
                    ei = ebest;
                    tcomp = tcomps(ie);
                end
                % Skip small error values for speed up
                if ids > 1 && sqrt(ei) < 0.02 * (ids-1)
                    continue
                end
                % Determine parameter index for given error value
                ipm = 1;
                for ipm1 = 1 : size(pln,4)
                    if pln(ids, ipc, ia, ipm1).emax <= ei
                        ipm = ipm1;
                    end
                end
                % GT n planes
                if fixngt
                    ipm = -1;
                    npln = -1;
                    ngt1 = max([sum((unique(pc{ids,ipc}.Intensity) >= 10)), 1]);
                    for ipm1 = 1 : size(pln,4)
                        plni1 = pln(ids, ipc, ia, ipm1).plane;
                        npln1 = sum(arrayfun(@(x) numel(x.index), plni1) > 1);
                        if npln1 >= ngt1
                            ipm = ipm1;
                            npln = npln1;
                        end
                    end
                    fprintf('ia = %i ngt = %i npln = %i ipm = %i\n', ia, ngt1, npln, ipm);
                    if ipm < 1
                        continue
                    end
                end
                % If same index as previously, just copy results and skip
                if ipmprev == ipm && ~scantcomp
                    nmsall(ids,ipc,ia,ie) =  nmsall(ids,ipc,ia,ie-1);
                    ncorrseg(ids,ipc,ia,ie) = ncorrseg(ids,ipc,ia,ie-1);
                    noverseg(ids,ipc,ia,ie) = noverseg(ids,ipc,ia,ie-1);
                    nundrseg(ids,ipc,ia,ie) = nundrseg(ids,ipc,ia,ie-1);
                    nmissseg(ids,ipc,ia,ie) = nmissseg(ids,ipc,ia,ie-1);
                    nnoisseg(ids,ipc,ia,ie) = nnoisseg(ids,ipc,ia,ie-1);
                    angdiff(ids,ipc,ia,ie) = angdiff(ids,ipc,ia,ie-1);
                    nangdiff(ids,ipc,ia,ie) = nangdiff(ids,ipc,ia,ie-1);
                    fcorrseg(ids,ipc,ia,ie) = fcorrseg(ids,ipc,ia,ie-1);
                    rmsems(ids,ipc,ia,ie) = rmsems(ids,ipc,ia,ie-1);
                    continue;
                end
                ipmprev = ipm;
                % Filter out atomic planes
                plni = rmatompln(pln(ids, ipc, ia, ipm));
                plni = plni.plane;
                npln = numel(plni);
                fprintf('%i ', npln);
                % SegComp evaluation
                [nmsall(ids,ipc,ia,ie), ...
                    ncorrseg(ids,ipc,ia,ie), ...
                    noverseg(ids,ipc,ia,ie), ...
                    nundrseg(ids,ipc,ia,ie), ...
                    nmissseg(ids,ipc,ia,ie), ...
                    nnoisseg(ids,ipc,ia,ie), ...
                    angdiff(ids,ipc,ia,ie), ...
                    nangdiff(ids,ipc,ia,ie), ...
                    fcorrseg(ids,ipc,ia,ie)] = ...
                    segcompeval(pci, plni, gtang{ids,ipc}, tcomp);
                % Calculate RMSE of reconstructed planes
                err = [];
                for ipln = 1 : npln
                    loc = locall(plni(ipln).index,:);
                    try
                        erri = calcerr(plni(ipln).param, loc, false);
                        rmsei = sqrt(mean(erri.^2));
                        if rmsei < 10
                            err(end+1:end+numel(erri)) = erri;
                        else
                            warning('Removing outlier plane with RMSE = %f for ids = %i, ia = %i, ipc = %i\n', rmsei, ids, ia, ipc);
                        end
                    catch
                        fprintf('Calc err failed\n');
                    end
                end
                rmsems(ids,ipc,ia,ie) = sqrt(mean(err.^2));
            end
            % Update progress display.
            it = it + 1;
            fprintf('Progress: %i/%i.\n', it, nit)            
        end
        % Calculate RMSE of ground-truth planes
        err = [];
        for igt = 1 : ngt
            % Select points
            loc = locall(pci.Intensity == gti(igt),:);
            % Initialize plane parameters based on extrema
            ex = minmax(loc');
            if ids == 1
                plninit = [
                    ex(1,1) ex(2,1) ex(3,2); ...
                    ex(1,2) ex(2,1) ex(3,2); ...
                    ex(1,1) ex(2,2) ex(3,2)];
            else
                plninit = [
                    ex(1,1) ex(2,1) ex(3,1); ...
                    ex(1,1) ex(2,1) ex(3,2); ...
                    ex(1,1) ex(2,2) ex(3,1)];
            end
            % Do the fit and calculate errors
            try
                [erri,xi] = calcerr(plninit, loc, true);
                err(end+1:end+numel(erri)) = erri;
            catch
                fprintf('Calc err failed\n');
            end
            if plotfits
                % Print and draw fit result
                fprintf('%i %i %f \n', gti(igt), numel(loc)/3, ...
                    sqrt(mean(erri.^2)));
                param = plninit./vecnorm(plninit,2,2).*xi;
                plane = struct('index', find(pci.Intensity == ...
                    gti(igt)),'error',[],'param',param);
                clf
                plotplane(pci, plane);
                campos([0 0.2 0]);
                drawnow;
            end
        end
        rmsegt(ids,ipc) = sqrt(mean(err.^2));
    end
end

%% Save evaluation.
evalfile = fullfile('result','eval.mat');
save(evalfile, 'ngtall', 'ncorrseg', 'noverseg', ...
    'nundrseg', 'nmissseg', 'nnoisseg', 'angdiff', 'nangdiff', ...
    'fcorrseg', 'nmsall', 'rmsems', 'rmsegt', 'algorithmname');
if scantcomp
    save(evalfile, 'tcomps', 'ebest', '-append')
else
    save(evalfile, 'tcomp', 'e', '-append')
end

end

function [err,x] = calcerr(pln, loc, dofit)
    xdata = pln;
    xdata(4:3+size(loc,1),:) = loc;
    xdata = double(xdata);
    rdata = vecnorm(xdata, 2, 2);
    xdata = xdata ./ rdata;
    x0 = rdata(1:3);
    if dofit
        cpuopt = optimoptions('lsqcurvefit', ...
            'Algorithm', 'levenberg-marquardt', ...
            'SpecifyObjectiveGradient', true, 'Display', 'none');
        ydata = rdata(4:end);
        [x,~] = lsqcurvefit(@rxp, x0, xdata, ydata, [], [], cpuopt);
    else
        x = x0;
    end
    [r,~] = rxp(x, xdata);
    err = rdata(4:end) - r;
end
