% Demo to fit multiple Gaussians to data.
clc;    % Clear the command window.
fprintf('Beginning to run %s.m.\n', mfilename);
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
clear global;
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

% First specify how many Gaussians there will be.
numGaussians = 6;

% Now make up some random parameters to create a set of Gaussians that we
% will sum together to get our test signal, from which we will try to guess the 
% parameters of the Gaussian curves that went into creating the test signal.
% Make centers in the range 0 - 100
centers = randi(100, 1, numGaussians);
% Make widths in the range 0 - 20
sigmas = randi(20, 1, numGaussians);
% Make amplitudes in the range 0 - 40
amplitudes = randi([10, 40], 1, numGaussians);
% Make signal that is the sum of all Gaussians
% g = gaussian(x, peakPosition, width)
x = linspace(0, 150, 1000);
y = zeros(1, length(x));
hFig = figure;

% Put all the parameters into a table for convenience in looking at, and using, the results.
tActual = table((1:numGaussians)', amplitudes(:), centers(:), sigmas(:), 'VariableNames', {'Number', 'Amplitude', 'Mean', 'Width'});
% Now sort parameters in order of increasing mean, just so it's easier to think about (though it's not required).
tActual = sortrows(tActual, 3);
tActual.Number = (1:numGaussians)' % Unsort the first column of numbers.

% Sum up the component curves to make our test signal that we will analyze to try to guess the component curves from.
legendStrings = cell(numGaussians, 1);
for k = 1 : numGaussians
	thisGaussian = tActual.Amplitude(k) * gaussian(x, tActual.Mean(k), tActual.Width(k));
	y = y + thisGaussian;
	plot(x, thisGaussian, '-', 'LineWidth', 1);
	hold on;
	legendStrings{k} = sprintf('Actual Gaussian %d', k);
	fprintf('Gaussian #%d has amplitude %5.1f, mean %5.1f, and sigma %5.1f.\n', k, tActual.Amplitude(k), tActual.Mean(k), tActual.Width(k));
end

% Optionasl: Add a tiny bit of noise.
noiseAmplitude = 0.03 * max(y);	% Add 3% noise.
y = y + noiseAmplitude * (rand(size(y)) - 0.5);

% Plot initial starting signal (the sum of the Gaussians).
hFig.WindowState = 'maximized';
hFig.Name = 'Original component curves summed together to form random test signal';
plot(x, y, 'k-', 'LineWidth', 2);
grid on
xlim(sort([x(1) x(end)]));
hold on
xlabel('x', 'FontSize', fontSize)
ylabel('y', 'FontSize', fontSize)
caption = sprintf('Sum of %d Gaussians', numGaussians);
title(caption, 'FontSize', fontSize, 'Interpreter', 'none');
legendStrings{end+1} = sprintf('Sum of all %d Gaussians', numGaussians);
legend(legendStrings);
drawnow;

%----------------------------------------------------------------------------------------------------------------------------------
% Now we have our test signal and we can begin....
% Fit Gaussian Peaks:
% Initial Gaussian Parameters
initialGuesses = [tActual.Mean(:), tActual.Width(:)];
% Add a little noise so that our first guess is not dead on accurate.
initialGuesses = initialGuesses + 2 * rand(size(initialGuesses))
startingGuesses = reshape(initialGuesses', 1, [])

global c NumTrials TrialError
% 	warning off

% Initializations
NumTrials = 0;  % Track trials
TrialError = 0; % Track errors
% t and y must be row vectors.
tFit = reshape(x, 1, []);
y = reshape(y, 1, []);

%-------------------------------------------------------------------------------------------------------------------------------------------
% Perform an iterative fit using the FMINSEARCH function to optimize the height, width and center of the multiple Gaussians.
options = optimset('TolX', 1e-4, 'MaxFunEvals', 10^12);  % Determines how close the model must fit the data
% First, set some options for fminsearch().
options.TolFun = 1e-4;
options.TolX = 1e-4;
options.MaxIter = 100000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEAVY LIFTING DONE RIGHT HERE:
% Run optimization
[parameter, fval, flag, output] = fminsearch(@(lambda)(fitgauss(lambda, tFit, y)), startingGuesses, options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%----------------------------------------------------------------------------------------------------------------
% Now plot results.
yhat = PlotComponentCurves(x, y, tFit, c, parameter);
% Compute the residuals between the actual y and the estimated y and put that into the graph's title.
meanResidual = mean(abs(y - yhat));
fprintf('The mean of the absolute value of the residuals is %f.\n', meanResidual);
caption = sprintf('Estimation of %d Gaussian Curves that will fit data.  Mean Residual = %f.', numGaussians, meanResidual);
title(caption, 'FontSize', fontSize, 'Interpreter', 'none');
drawnow;

% Make table for the fitted, estimated results.
% First make numGaussians row by 3 column matrix: Column 1 = amplitude, column 2 = mean, column 3 = width.
% 	parameter % Print to command window.
estimatedMuSigma = reshape(parameter, 2, [])';
gaussianParameters = [c, estimatedMuSigma];
% Now sort parameters in order of increasing mean
gaussianParameters = sortrows(gaussianParameters, 2);
tActual % Display actual table in the command window.
% Create table of the output parameters and display it below the actual, true parameters.
tEstimate = table((1:numGaussians)', c(:), estimatedMuSigma(:, 1), estimatedMuSigma(:, 2), 'VariableNames', {'Number', 'Amplitude', 'Mean', 'Width'})

% Plot the error as a function of trial number.
hFigError = figure();
hFigError.Name = 'Errors';
plot(TrialError, 'b-');
% hFigError.WindowState = 'maximized';
grid on;
xlabel('Trial Number', 'FontSize', fontSize)
ylabel('Error', 'FontSize', fontSize)

caption = sprintf('Errors for all %d trials.', length(TrialError));
title(caption, 'FontSize', fontSize, 'Interpreter', 'none');

message = sprintf('Done!\nHere is the result!\nNote: there could be multiple ways\n(mulitple sets of Gaussians)\nthat you could achieve the same sum (same test curve).');
fprintf('Done running %s.m.\n', mfilename);
msgboxw(message);

%=======================================================================================================================================================
function yhat = PlotComponentCurves(x, y, t, c, parameter)
try
	fontSize = 20;
	% Get the means and widths.
	means = parameter(1 : 2 : end);
	widths = parameter(2 : 2 : end);
	% Now plot results.
	hFig2 = figure;
	hFig2.Name = 'Fitted Component Curves';
% 	plot(x, y, '--', 'LineWidth', 2)
	hold on;
	yhat = zeros(1, length(t));
	numGaussians = length(c);
	legendStrings = cell(numGaussians + 2, 1);
	for k = 1 : numGaussians
		% Get each component curve.
		thisEstimatedCurve = c(k) .* gaussian(t, means(k), widths(k));
		% Plot component curves.
		plot(x, thisEstimatedCurve, '-', 'LineWidth', 2);
		hold on;
		% Overall curve estimate is the sum of the component curves.
		yhat = yhat + thisEstimatedCurve;
		legendStrings{k} = sprintf('Estimated Gaussian %d', k);
	end
	% Plot original summation curve, that is the actual curve.
	plot(x, y, 'r-', 'LineWidth', 1)
	% Plot estimated summation curve, that is the estimate of the curve.
	plot(x, yhat, 'k--', 'LineWidth', 2)
	grid on;
	xlabel('X', 'FontSize', fontSize)
	ylabel('Y', 'FontSize', fontSize)
	caption = sprintf('Estimation of %d Gaussian Curves that will fit data.', numGaussians);
	title(caption, 'FontSize', fontSize, 'Interpreter', 'none');
	grid on
	legendStrings{numGaussians+1} = sprintf('Actual original signal');
	legendStrings{numGaussians+2} = sprintf('Sum of all %d Gaussians', numGaussians);
	legend(legendStrings);
	xlim(sort([x(1) x(end)]));
	hFig2.WindowState = 'maximized';
	drawnow;
	
catch ME
	% Some error happened if you get here.
	callStackString = GetCallStack(ME);
	errorMessage = sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', ...
		mfilename, callStackString, ME.message);
	WarnUser(errorMessage);
end
end % of PlotComponentCurves


%=======================================================================================================================================================
function theError = fitgauss(lambda, t, y)
% Fitting function for multiple overlapping Gaussians, with statements
% added (lines 18 and 19) to slow the progress and plot each step along the
% way, for educational purposes.
% Author: T. C. O'Haver, 2006

global c NumTrials TrialError
try
	
	A = zeros(length(t), round(length(lambda) / 2));
	for j = 1 : length(lambda) / 2
		A(:,j) = gaussian(t, lambda(2 * j - 1), lambda(2 * j))';
	end
	
	c = A \ y';
	z = A * c;
	theError = norm(z - y');
	
	% Penalty so that heights don't become negative.
	if sum(c < 0) > 0
		theError = theError + 1000000;
	end
	
	NumTrials = NumTrials + 1;
	TrialError(NumTrials) = theError;
catch ME
	% Some error happened if you get here.
	callStackString = GetCallStack(ME);
	errorMessage = sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', ...
		mfilename, callStackString, ME.message);
	WarnUser(errorMessage);
end
end % of fitgauss()


%=======================================================================================================================================================
function g = gaussian(x, peakPosition, width)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x - peakPosition) ./ (0.60056120439323 .* width)) .^ 2);
end % of gaussian()