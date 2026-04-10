function runDynareWithCache(modFile, macroArgs, extraArgs)
% runDynareWithCache  Run Dynare, skipping preprocessing if cached output exists.
%
%   runDynareWithCache(modFile, macroArgs, extraArgs)
%
%   modFile   - Name of the .mod file (e.g., 'dynamicModel.mod')
%   macroArgs - String of -D flags (e.g., '-DnMeasure=2 -Dperturbation_order=1')
%   extraArgs - String of other Dynare options (e.g., 'noclearall nograph')
%               (optional, default '')
%
% Caching mechanism:
%   After the first Dynare call, the preprocessor output is saved in
%   +MODELNAME/. Subsequent calls with IDENTICAL arguments (both macro
%   and extra) will skip preprocessing and call MODELNAME.driver directly.
%   If arguments change, a full dynare call is performed.
%
% The cache key is stored in MODELNAME_cache_key.mat alongside the
% preprocessed output.

if nargin < 3
    extraArgs = '';
end

% Extract model name (strip .mod/.dyn extension)
[~, modelName] = fileparts(modFile);

% Build the full argument string as the cache key.
% All args (macro and extra) affect the preprocessor output, so include all.
cacheKey = strtrim([extraArgs ' ' macroArgs]);

% Build the cache key file path
cacheKeyFile = [modelName '_cache_key.mat'];

% Check if preprocessed output exists
driverExists = exist(['+' modelName filesep 'driver.m'], 'file') == 2;

% Check if cache key matches
cacheHit = false;
if driverExists && exist(cacheKeyFile, 'file') == 2
    cached = load(cacheKeyFile, 'cachedKey');
    if isfield(cached, 'cachedKey') && strcmp(cached.cachedKey, cacheKey)
        cacheHit = true;
    end
end

if cacheHit
    % Cache hit: skip preprocessing, run driver directly
    fprintf('Dynare cache hit for %s [%s]\n', modFile, cacheKey);
    fprintf('Skipping preprocessing, running driver directly...\n');

    % Ensure Dynare is configured (adds paths, etc.)
    dynare_config();

    % Clear the driver from MATLAB's cache to pick up any file changes
    clear(['+' modelName '/driver']);

    % Run the driver in the caller's workspace (same as dynare.m line 306)
    evalin('caller', [modelName '.driver']);
else
    % Cache miss: run full dynare command
    fprintf('Dynare cache miss for %s [%s] — running full preprocessing...\n', modFile, cacheKey);

    % Build the full dynare command
    cmd = ['dynare ' modFile ' ' extraArgs ' ' macroArgs];
    evalin('caller', cmd);

    % Save cache key for future runs
    cachedKey = cacheKey;
    save(cacheKeyFile, 'cachedKey');
end

end
