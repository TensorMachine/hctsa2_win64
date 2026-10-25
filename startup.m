function startup()
% STARTUP   Add all paths required for the hctsa package.

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% Resolve paths relative to this file's own location (not the current
% working directory), so startup works regardless of where it's called from:
hctsaDir = fileparts(mfilename('fullpath'));

% We use this function a bit:
addfcn = @(x) addpath(fullfile(hctsaDir,x));

catch22Dir = fullfile(hctsaDir,'Toolboxes','catch22','wrap_Matlab');
% hctsa2: catch22 is vendored in Toolboxes/catch22 (not a git submodule), so
% upstream's automatic "git submodule update" is not needed and is left out: it
% would also need git and a network connection.

fprintf(1,'Adding paths for the highly comparative time-series analysis package...\n')

% ------------------------------------------------------------------------------
%% First add all the basic paths:
% ------------------------------------------------------------------------------
% hctsa2: the root itself, which holds the HCTSA2_* tools (HCTSA2_Analyse,
% HCTSA2_CheckInstall, HCTSA2_BuildTisean, ...). They resolve from any folder
% once startup has run, e.g. when it is called from a MATLAB startup file.
addpath(hctsaDir);
addfcn('FeatureSets'); % files for setting input files corresponding to various feature sets
addfcn('Calculation'); % code for calculating results
addfcn('PlottingAnalysis'); % code for analysing and plotting results
addfcn('Operations'); % core code files for performing operations
addfcn('PeripheryFunctions'); % periphery functions used in the code toolbox
addfcn('TimeSeries'); % time series data files for analysis
% hctsa2 (3.0.0): Tests/ holds HCTSA2's test and debug tools (HCTSA2_TestRun,
% HCTSA2_TestTisean, ...) as well as upstream's unit tests, so the top level holds
% only the tools a user works with; on the path, the test tools run by name as before.
addfcn('Tests');
fprintf(1,'Core directories added.\n')

% ------------------------------------------------------------------------------
%% Now add all the external code packages and periphery toolboxes:
% ------------------------------------------------------------------------------
fprintf(1,'Adding external time-series toolboxes...')
% Kaplan's routines:
fprintf(1,' Danny Kaplan')
addpath(fullfile(hctsaDir,'Toolboxes','Danny_Kaplan'));

% Gaussian Process Toolbox, gpml, by Carl Edward Rasmussen & Hannes Nickisch:
fprintf(1,', Gaussian Process Code\n')
addpath(fullfile(hctsaDir,'Toolboxes','gpml'));
GP_startup % add nested directories

% Zoubin Gharamani's hmm toolbox, ZG_hmm
fprintf(1,'HMM toolbox')
addpath(fullfile(hctsaDir,'Toolboxes','ZG_hmm'));

% ARFIT Toolbox
fprintf(1,', ARfit toolbox')
addpath(fullfile(hctsaDir,'Toolboxes','ARFIT'));

% Michael Small's utilities
fprintf(1,', Michael Small\n')
addpath(fullfile(hctsaDir,'Toolboxes','Michael_Small'));

% DVV Toolbox
fprintf(1,'DVV Toolbox')
addpath(fullfile(hctsaDir,'Toolboxes','DVV_Toolbox'));

% Physionet
fprintf(1,', Physionet');
addpath(fullfile(hctsaDir,'Toolboxes','Physionet'));

% Max Little's steps/bumps toolbox
fprintf(1,', Max Little''s steps_bumps toolkit')
addpath(fullfile(hctsaDir,'Toolboxes','Max_Little','steps_bumps_toolkit'));

% Max Little's rpde code
fprintf(1,', rpde')
addpath(fullfile(hctsaDir,'Toolboxes','Max_Little','rpde'));

% nsamdf
fprintf(1,', nsamdf,\n');
addpath(fullfile(hctsaDir,'Toolboxes','nsamdf'));

% catch22 (vendored; optional: no operation of the hctsa library needs it)
if isfolder(catch22Dir)
    fprintf(1,'catch22, ')
    addpath(catch22Dir);
end

% Kraskov (KSG) mutual information, the KSG_mi MEX (hctsa2 2.4.0: replaces the Java
% information dynamics toolkit, JIDT, which hctsa2 no longer ships)
fprintf(1,'KSG mutual information')
addpath(fullfile(hctsaDir,'Toolboxes','KSG'));

% ------------------------------------------------------------------------------
% Add paths for the TISEAN and ripser executables.
% ------------------------------------------------------------------------------
% hctsa2: platform-aware. On Windows the executables are built by
% HCTSA2_BuildTisean into Toolboxes/Tisean_bin (ripser.exe included) and PATH
% entries are separated by ';'. Elsewhere, upstream's layout: TISEAN compiled
% into Toolboxes/Tisean_3.0.1/bin by compile_tisean (falling back to ~/bin for
% older manual installations), ripser into Toolboxes/ripser/bin. Prepended, not
% appended, so an older install elsewhere on the PATH cannot shadow them.
% BF_RunTisean looks in the same places directly, so this is for the shell.
fprintf(1,', TISEAN')
if ispc
    tiseanBinaryLocation = fullfile(hctsaDir,'Toolboxes','Tisean_bin');
else
    tiseanBinaryLocation = fullfile(hctsaDir,'Toolboxes','Tisean_3.0.1','bin');
    if ~isfolder(tiseanBinaryLocation)
        homeDir = getenv('HOME');
        if isempty(homeDir)
            [~,homeDir] = system('echo $HOME'); % get system home directory
            homeDir = regexprep(homeDir,'[\s]',''); % remove whitespace
        end
        tiseanBinaryLocation = fullfile(homeDir,'bin');
    end
end
ripserBinaryLocation = fullfile(hctsaDir,'Toolboxes','ripser','bin');
fprintf(1,', ripser\n')
for binDir = {ripserBinaryLocation,tiseanBinaryLocation}
    if isfolder(binDir{1}) && ~any(strcmp(strsplit(getenv('PATH'),pathsep),binDir{1}))
        setenv('PATH',[binDir{1},pathsep,getenv('PATH')]);
    end
end
if isfolder(tiseanBinaryLocation)
    fprintf(1,'TISEAN binaries: %s\n',tiseanBinaryLocation);
else
    fprintf(1,['TISEAN binaries not found (%s) -- TISEAN-based features will error. ' ...
        'Build them with HCTSA2_BuildTisean (Windows) or compile_tisean.\n'],tiseanBinaryLocation);
end

% ------------------------------------------------------------------------------
%% Finished:
% ------------------------------------------------------------------------------
fprintf(1,'\n---Done.\n')

end
