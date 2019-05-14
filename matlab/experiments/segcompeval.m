function [nmsall, ncorrseg, noverseg, nundrseg, nmissseg, nnoisseg, ...
    angdiff, nangdiff, fcorrseg, corrseg] = ...
    segcompeval(pc, pln, gtang, tcomp)
% SEGCOMPEVAL Evaluate plane reconstruction performance.
%   Calculates metrics from the SegComp paper:
%   "An Experimental Comparison of Range Image Segmentation Algorithms,"
%   in IEEE Transactions on Pattern Analysis and Machine Intelligence,
%   July 1996, pp. 1-17

% Copyright Daniel Buescher 2018

%% Calculate overlap matrix
npln = numel(pln);
gti = unique(pc.Intensity);
% Only ground truth indices >= 10 are valid:
gti = gti(gti >= 10);
ngt = numel(gti);
pgt = NaN(ngt, 1);
pms = NaN(npln, 1);
o = NaN(ngt, npln);
for igt = 1 : ngt % number of gt-planes
    for ipln = 1 : npln % number of extracted planes
        indpln = pln(ipln).index; % p-indices of extracted planes
        indgt = find(pc.Intensity == gti(igt)); % indices of corresponding gt plane
        pms(ipln, 1) = numel(indpln); % number of extracted points
        pgt(igt, 1) = numel(indgt); % number of gt points
        o(igt, ipln) = numel(intersect(indpln, indgt)); % number of intersected points
    end
end

%% Calculate fractions of overlapping points
fms = o ./ pms';
fgt = o ./ pgt;
foverseg = sum((fms>tcomp) .* fgt, 2) .* (sum(fms>tcomp, 2) >= 2);
fundrseg = sum((fgt>tcomp) .* fms, 1) .* (sum(fgt>tcomp, 1) >= 2);
% Calculate segmentation types
corrseg = (fms > tcomp) & (fgt > tcomp);
overseg = foverseg > tcomp;
undrseg = fundrseg > tcomp;
% Calculate average fractions
avgcorrseg = (fms .* corrseg + fgt .* corrseg) / 2;
avgoverseg = (sum(fms .* (fms > tcomp), 2) + foverseg .* overseg) ...
    ./ (sum(fms > tcomp, 2) + 1) .* overseg;
avgundrseg = (sum(fgt .* (fgt > tcomp), 1) + fundrseg .* undrseg) ...
    ./ (sum(fgt > tcomp, 1) + 1) .* undrseg;
avgoversegcause = sum((fms>tcomp) .* overseg .* fms, 1);
avgundrsegcause = sum((fgt>tcomp) .* undrseg .* fgt, 2);
% Exclusify by taking larger average fractions
corrseg = (avgcorrseg > avgoverseg & avgcorrseg > avgundrseg);
overseg = avgoverseg > sum(avgcorrseg, 2) & avgoverseg > avgundrsegcause;
undrseg = avgundrseg > sum(avgcorrseg, 1) & avgundrseg > avgoversegcause;
% Find remaining planes
oversegcause = sum((fms>tcomp) .* overseg, 1);
undrsegcause = sum((fgt>tcomp) .* undrseg, 2);
missseg = ~(sum(corrseg, 2) + overseg + undrsegcause);
noisseg = ~(sum(corrseg, 1) + undrseg + oversegcause);
% Calculate metrics
nmsall = npln;
ncorrseg = sum(sum(corrseg));
noverseg = sum(overseg);
nundrseg = sum(undrseg);
nmissseg = sum(missseg);
nnoisseg = sum(noisseg);

% Consistency check
testgt = [sum(corrseg, 2), overseg, undrsegcause, missseg];
testms = [sum(corrseg, 1); undrseg; oversegcause; noisseg];
if ~isequal(sum(testgt, 2), ones(ngt, 1))
    warning('Ground-truth classification is not consistent!');
end
if ~isequal(sum(testms, 1), ones(1, npln))
    warning('Measurement classification is not consistent!');
end

% Loop over pairs of ground-thruth planes and compare the
% angle between them to the reconstucted ones (if a
% correct correpsondence exists)
angdiffsum = 0;
nangdiff1 = 0;
for iang = 1:size(gtang, 1)
    % Find ground thruth and reconstructed indices
    igt1 = gtang(iang,1);
    igt2 = gtang(iang,2);
    ims1 = find(corrseg(igt1-9,:));
    ims2 = find(corrseg(igt2-9,:));
    if ~isempty(ims1) && ~isempty(ims2)
        % Calculate reconstructed plane normals
        n1 = cross( ...
            pln(ims1).param(1,:) - pln(ims1).param(3,:), ...
            pln(ims1).param(2,:) - pln(ims1).param(3,:));
        n2 = cross( ...
            pln(ims2).param(1,:) - pln(ims2).param(3,:), ...
            pln(ims2).param(2,:) - pln(ims2).param(3,:));
        % Calculate angles
        ang = atan2(norm(cross(n1,n2)),dot(n1,n2))/pi*180;
        gtangi = gtang(iang,3);
        if ang > 90
            ang = 180 - ang;
        end
        if gtangi > 90
            gtangi = 180 - gtangi;
        end
        angdiffsum = angdiffsum + abs(ang-gtangi);
        nangdiff1 = nangdiff1 + 1;
    end
end
angdiff = angdiffsum;
nangdiff = nangdiff1;

gtcorr = sum(corrseg, 2);
ngti = arrayfun(@(x) sum(pc.Intensity(:) == x), gti);
fcorrseg = sum(ngti.*gtcorr) / sum(ngti);

end
