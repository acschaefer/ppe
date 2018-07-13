function planetorays()
% Example of the Matlab binding of the Gpufit library implementing
% Levenberg Marquardt curve fitting in CUDA
% https://github.com/gpufit/Gpufit
%
% Multiple fits of a 2D Gaussian peak function with Poisson distributed noise
% http://gpufit.readthedocs.io/en/latest/bindings.html#matlab

assert(gpufit_cuda_available(), 'CUDA not available');

%% number of fits and fit points
number_fits = 1e4;
size_x = 20;
number_parameters = 2;

%% set input arguments

% true parameters
true_parameters = single([10, 5]);

% initialize random number generator
rng(0);

% initial parameters (can be randomized)
initial_parameters = repmat(single(ones(size(true_parameters'))), [1, number_fits]);

% generate x values
g = single(0 : size_x - 1);
x = ndgrid(g);

% generate data with Poisson noise
data = linear_1d(x, true_parameters);
data = repmat(data(:), [1, number_fits]);
data = poissrnd(data);

% tolerance
tolerance = 1e-3;

% maximum number of iterations
max_n_iterations = 20;

% estimator id
estimator_id = EstimatorID.MLE;

% model ID
model_id = ModelID.PLANETORAYS;
%model_id = ModelID.LINEAR_1D;

%% run Gpufit
[parameters, states, chi_squares, n_iterations, time] = gpufit(data, [], ...
    model_id, initial_parameters, tolerance, max_n_iterations, [], estimator_id, []);

%% displaying results
display_results('Plane to rays', model_id, number_fits, number_parameters, size_x, time, true_parameters, parameters, states, chi_squares, n_iterations);

end


function g = linear_1d(x, p)
g = p(1) + p(2) * x;
end

function display_results(name, model_id, number_fits, number_parameters, size_x, time, true_parameters, parameters, states, chi_squares, n_iterations)

%% displaying results
converged = states == 0;
fprintf('\nGpufit of %s\n', name);

% print summary
fprintf('\nmodel ID:        %d\n', model_id);
fprintf('number of fits:  %d\n', number_fits);
fprintf('fit size:        %d x %d\n', size_x, size_x);
fprintf('mean chi-square: %6.2f\n', mean(chi_squares(converged)));
fprintf('mean iterations: %6.2f\n', mean(n_iterations(converged)));
fprintf('time:            %6.2f s\n', time);

% get fit states
number_converged = sum(converged);
fprintf('\nratio converged         %6.2f %%\n', number_converged / number_fits * 100);
fprintf('ratio max it. exceeded  %6.2f %%\n', sum(states == 1) / number_fits * 100);
fprintf('ratio singular hessian  %6.2f %%\n', sum(states == 2) / number_fits * 100);
fprintf('ratio neg curvature MLE %6.2f %%\n', sum(states == 3) / number_fits * 100);

% mean and std of fitted parameters
converged_parameters = parameters(:, converged);
converged_parameters_mean = mean(converged_parameters, 2);
converged_parameters_std  = std(converged_parameters, [], 2);
fprintf('\nparameters of %s\n', name);
for i = 1 : number_parameters
    fprintf('p%d true %6.2f mean %6.2f std %6.2f\n', i, true_parameters(i), converged_parameters_mean(i), converged_parameters_std(i));
end

end