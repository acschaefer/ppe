function planetorays()
% Example of the Matlab binding of the Gpufit library implementing
% Levenberg Marquardt curve fitting in CUDA
% https://github.com/gpufit/Gpufit
%
% Multiple fits of a 2D Gaussian peak function with Poisson distributed noise
% http://gpufit.readthedocs.io/en/latest/bindings.html#matlab

assert(gpufit_cuda_available(), 'CUDA not available');

%% number of fits and fit points
number_fits = 1e5;
number_parameters = 3;

%% set input arguments

% true parameters
true_parameters = double([10, 20, 20]);

% initialize random number generator
rng(0);

% initial parameters
initial_parameters = repmat(double(ones(size(true_parameters'))), [1, number_fits]);

% generate data
data = double([10; 20; 30; 10]);
data = repmat(data(:), [1, number_fits]);

% grow data linearly for validation purposes
data = data .* (1 + 0.001 * linspace(1, number_fits, number_fits));

% add Poisson noise to data
%data = poissrnd(data);

% 3 plane vectors and 3+1 ray vectors (fixed direction)
%user_info = double([1 0 0 0 1 0 0 0 1 1 1 1]);
%user_info = double([1 0 0 0 1 0 0 0 1 1 0 0]);
%user_info = repmat([user_info(1:9), user_info], [1, number_fits]);

% vectors with randomized directions
for i = 1:number_fits
   r1 = rand();
   r2 = rand();
   r3 = rand();
   r1a = sqrt(1 - r1*r1);
   r2a = sqrt(1 - r2*r2);
   r3a = sqrt(1 - r3*r3);
   vecs = double([r1 r1a 0 0 r2 r2a r3a 0 r3 r1 r1a 0]);
   user_info(:,i) = [vecs(1:9), vecs]';
end

% tolerance
tolerance = double(1e-6);

% maximum number of iterations
max_n_iterations = 20;

% estimator id
estimator_id = EstimatorID.LSE;

% model ID
model_id = ModelID.PLANETORAYS;

%% run Gpufit
[parameters, states, chi_squares, n_iterations, time] = gpufit(data, [], ...
    model_id, initial_parameters, tolerance, max_n_iterations, [], estimator_id, user_info);

%% displaying results
display_results('Plane to rays', model_id, number_fits, number_parameters, time, true_parameters, parameters, states, chi_squares, n_iterations);

end

function display_results(name, model_id, number_fits, number_parameters, time, true_parameters, parameters, states, chi_squares, n_iterations)

%% displaying results
converged = states == 0;
fprintf('\nGpufit of %s\n', name);

% print summary
fprintf('\nmodel ID:        %d\n', model_id);
fprintf('number of fits:  %d\n', number_fits);
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