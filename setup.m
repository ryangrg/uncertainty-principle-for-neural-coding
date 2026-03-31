% setup.m - Add all project directories to the MATLAB path
%
% Run this script before executing any simulations:
%   >> run('setup.m')

projectRoot = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(projectRoot, 'lib')));
addpath(fullfile(projectRoot, 'data'));

fprintf('Uncertainty Principle project paths added.\n');
fprintf('Project root: %s\n', projectRoot);
