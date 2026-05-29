function result = bootstrapROCByPatient(scores, labels, patientID, nBoot)
% BOOTSTRAPROCBYPATIENT
% Patient-level bootstrap ROC analysis for clustered data.
%
% Inputs
%   scores    : Nx1 numeric vector
%   labels    : Nx1 binary vector (0/1)
%   patientID : Nx1 vector (numeric or cellstr/string/categorical)
%   nBoot     : number of bootstrap replicates
%
% Output
%   result : struct with fields
%       .AUC              - apparent/full-sample AUC
%       .AUC_boot         - bootstrap AUCs
%       .AUC_mean         - mean bootstrap AUC
%       .AUC_CI           - 95% percentile CI
%       .FPR              - full-sample ROC x
%       .TPR              - full-sample ROC y
%       .Thresholds       - thresholds from perfcurve
%       .nPatients        - number of unique patients
%       .nSamples         - number of total observations

if nargin < 4 || isempty(nBoot)
    nBoot = 1000;
end

% force column vectors
scores = scores(:);
labels = labels(:);
patientID = patientID(:);

% basic checks
if numel(scores) ~= numel(labels) || numel(scores) ~= numel(patientID)
    error('scores, labels, and patientID must have the same length.');
end

validMask = ~(isnan(scores) | isnan(labels));
scores = scores(validMask);
labels = labels(validMask);
patientID = patientID(validMask);

% labels should be binary
uLabels = unique(labels);
if ~all(ismember(uLabels, [0 1]))
    error('labels must be binary and coded as 0/1.');
end

% full-sample ROC
[FPR, TPR, Thresholds, AUC] = perfcurve(labels, scores, 1);

% unique patients
[uniquePatients, ~, patientIndex] = unique(patientID, 'stable');
nPatients = numel(uniquePatients);

AUC_boot = nan(nBoot, 1);

for b = 1:nBoot
    % resample patient indices with replacement
    sampledPatientIdx = randi(nPatients, nPatients, 1);

    % build bootstrap sample by concatenating all segments
    bootScores = [];
    bootLabels = [];

    for j = 1:nPatients
        thisIdx = sampledPatientIdx(j);
        obsMask = (patientIndex == thisIdx);

        bootScores = [bootScores; scores(obsMask)];
        bootLabels = [bootLabels; labels(obsMask)];
    end

    % skip degenerate samples with only one class
    if numel(unique(bootLabels)) < 2
        continue;
    end

    [~, ~, ~, AUC_boot(b)] = perfcurve(bootLabels, bootScores, 1);
end

% remove invalid replicates
AUC_boot = AUC_boot(~isnan(AUC_boot));

result = struct();
result.AUC = AUC;
result.AUC_boot = AUC_boot;
result.AUC_mean = mean(AUC_boot);
result.AUC_CI = prctile(AUC_boot, [2.5 97.5]);
result.FPR = FPR;
result.TPR = TPR;
result.Thresholds = Thresholds;
result.nPatients = nPatients;
result.nSamples = numel(scores);
end