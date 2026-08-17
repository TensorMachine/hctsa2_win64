function startup()
% STARTUP   Add all paths required for the hctsa package.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

% We use this function a bit:
addfcn = @(x) addpath(fullfile(pwd,x));

fprintf(1,'Adding paths for the highly comparative time-series analysis package...\n')

% ------------------------------------------------------------------------------
%% First add all the basic paths:
% ------------------------------------------------------------------------------
addfcn('Database'); % code for setting up and communicating with the mySQL database
addfcn('FeatureSets'); % files for setting input files corresponding to various feature sets
addfcn('Calculation'); % code for calculating results
addfcn('PlottingAnalysis'); % code for analysing and plotting results
addfcn('Operations'); % core code files for performing operations
addfcn('PeripheryFunctions'); % periphery functions used in the code toolbox
addfcn('TimeSeries'); % time series data files for analysis
fprintf(1,'Core directories added.\n')

% ------------------------------------------------------------------------------
%% Now add all the external code packages and periphery toolboxes:
% ------------------------------------------------------------------------------
fprintf(1,'Adding external time-series toolboxes...')
% Kaplan's routines:
fprintf(1,' Danny Kaplan')
addpath(fullfile(pwd,'Toolboxes','Danny_Kaplan'));

% Code by Marwan, from CRP Toolbox version 5.17  (R28.16)
fprintf(1,', Marwan')
addpath(fullfile(pwd,'Toolboxes','Marwan_crptool'));

% Gaussian Process Toolbox, gpml, by Carl Edward Rasmussen & Hannes Nickisch:
fprintf(1,', Gaussian Process Code\n')
addpath(fullfile(pwd,'Toolboxes','gpml'));
GP_startup % add nested directories

% Zoubin Gharamani's hmm toolbox, ZG_hmm
fprintf(1,'HMM toolbox')
addpath(fullfile(pwd,'Toolboxes','ZG_hmm'));

% ARFIT Toolbox
fprintf(1,', ARfit toolbox')
addpath(fullfile(pwd,'Toolboxes','ARFIT'));

% Michael Small's utilities
fprintf(1,', Michael Small')
addpath(fullfile(pwd,'Toolboxes','Michael_Small'));

% Code from Matlab Central
fprintf(1,', Matlab Central code')
addpath(fullfile(pwd,'Toolboxes','MatlabCentral'));

% Rudy Moddemeijer's code
fprintf(1,', Rudy Moddemeijer\n')
addpath(fullfile(pwd,'Toolboxes','Rudy_Moddemeijer'));

% DVV Toolbox
fprintf(1,'DVV Toolbox')
addpath(fullfile(pwd,'Toolboxes','DVV_Toolbox'));

% Physionet
fprintf(1,', Physionet');
addpath(fullfile(pwd,'Toolboxes','Physionet'));

% Max Little's steps/bumps toolbox
fprintf(1,', Max Little''s steps_bumps toolkit')
addpath(fullfile(pwd,'Toolboxes','Max_Little','steps_bumps_toolkit'));

% Max Little's fastdfa code
fprintf(1,', fastdfa')
addpath(fullfile(pwd,'Toolboxes','Max_Little','fastdfa'));

% Max Little's rpde code
fprintf(1,', rpde')
addpath(fullfile(pwd,'Toolboxes','Max_Little','rpde'));

% nsamdf
fprintf(1,', nsamdf,\n');
addpath(fullfile(pwd,'Toolboxes','nsamdf'));

% Misc code
fprintf(1,'misc')
addpath(fullfile(pwd,'Toolboxes','Misc'));

% catch22
fprintf(1,', catch22')
addpath(fullfile(pwd,'Toolboxes','catch22','wrap_Matlab'));

% TSTOOL was removed in upstream v2.0.0 and its operations replaced with native
% and TISEAN-based equivalents. Nothing adds OpenTSTOOL to the path any more --
% the toolbox itself has been deleted from Toolboxes/.

% Java information dynamics toolkit written by Joseph Lizier
% (should be ok to re-add this every time startup is run)
fprintf(1,', Information dynamics toolkit, ')
javaaddpath(fullfile(pwd,'Toolboxes','infodynamics-dist','infodynamics.jar'));

% ------------------------------------------------------------------------------
% Add path for TISEAN binaries.
% ------------------------------------------------------------------------------
% Platform-aware. The original assumed macOS/Linux: it read $HOME via a shell
% echo, joined PATH with ':', and set DYLD_LIBRARY_PATH (a macOS-only variable).
% All three are wrong on Windows, where the separator is ';' and the TISEAN
% executables are .exe files that hctsa 2.0 expects under Toolboxes/Tisean_bin.
if ispc
    tiseanBinaryLocation = fullfile(pwd,'Toolboxes','Tisean_bin');
    pathSep = ';';
else
    homeDir = getenv('HOME');
    if isempty(homeDir)
        [~,homeDir] = system('echo $HOME');
        homeDir = regexprep(homeDir,'[\s]','');
    end
    tiseanBinaryLocation = fullfile(homeDir,'bin');
    pathSep = ':';
end
if exist(tiseanBinaryLocation,'dir') == 7
    if isempty(strfind(getenv('PATH'),tiseanBinaryLocation)) %#ok<STREMP>
        setenv('PATH',[getenv('PATH'),pathSep,tiseanBinaryLocation]);
    end
    fprintf(1,'TISEAN binaries: %s\n',tiseanBinaryLocation);
else
    fprintf(1,'TISEAN binaries not found (%s) -- TISEAN-based features will error.\n', ...
        tiseanBinaryLocation);
end

% DYLD_LIBRARY_PATH is macOS-specific; setting it elsewhere is harmless but
% pointless, and on Windows it is meaningless.
if ismac
    setenv('DYLD_LIBRARY_PATH','/usr/local/bin');
end

% ------------------------------------------------------------------------------
%% Finished:
% ------------------------------------------------------------------------------
fprintf(1,'\n---Done.\n')

end
