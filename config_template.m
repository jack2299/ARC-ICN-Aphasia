%% CONFIGURATION TEMPLATE for ARC Pipeline
% =========================================================================
% INSTRUCTIONS:
% 1. Copy this file to 'config_local.m' in the same directory
% 2. Edit 'config_local.m' with your actual paths
% 3. DO NOT commit 'config_local.m' to GitHub (it's in .gitignore)
% =========================================================================

function config = config_template()
    % =====================================================================
    % REQUIRED: Root directory containing your data
    % =====================================================================
    % Example: config.rootDir = '/path/to/your/processed/files';
    config.rootDir = '';  % <<< EDIT THIS IN config_local.m
    
    % =====================================================================
    % OPTIONAL: Override individual paths (if different from defaults)
    % Defaults are derived from rootDir unless specified here
    % =====================================================================
    config.spmDir = [];        % SPM.mat files location (default: fullfile(rootDir, 'SPM_mat_files_backup'))
    config.clinicalFile = [];  % Clinical data file (default: fullfile(rootDir, 'participants.txt'))
    config.outputDir = [];     % Output directory (default: fullfile(rootDir, 'Analysis_Output'))
    
    % =====================================================================
    % PIPELINE SETTINGS (optional - override defaults if needed)
    % =====================================================================
    config.verbose = true;     % Print progress messages
    config.saveCSV = true;     % Save CSV outputs
    config.saveMAT = true;     % Save MAT outputs
    
    % =====================================================================
    % VALIDATION: Check if required fields are set
    % =====================================================================
    if isempty(config.rootDir)
        error('config_local.m: rootDir must be set. Edit config_local.m with your data path.');
    end
    
    % Set default paths if not overridden
    if isempty(config.spmDir)
        config.spmDir = fullfile(config.rootDir, 'SPM_mat_files_backup');
    end
    if isempty(config.clinicalFile)
        config.clinicalFile = fullfile(config.rootDir, 'participants.txt');
    end
    if isempty(config.outputDir)
        config.outputDir = fullfile(config.rootDir, 'Analysis_Output');
    end
end
