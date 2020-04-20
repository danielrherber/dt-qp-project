%--------------------------------------------------------------------------
% INSTALL_DTQP
% Project link: https://github.com/danielrherber/dt-qp-project
% This scripts helps you get the DT QP Project up and running
%--------------------------------------------------------------------------
% Automatically adds project files to your MATLAB path, downloads the
% required files, and opens an example
%--------------------------------------------------------------------------
% Install script based on MFX Submission Install Utilities
% https://github.com/danielrherber/mfx-submission-install-utilities
% https://www.mathworks.com/matlabcentral/fileexchange/62651
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber)
% Additional contributors:
% - Yong Hoon Lee (yonghoonlee)
% - Athul K. Sundarrajan (AthulKrishnaSundarrajan)
%--------------------------------------------------------------------------
function INSTALL_DTQP(varargin)

% initialize
silentflag = 0; % don't be silent

% parse inputs
if ~isempty(varargin)
    if any(strcmpi(varargin,'silent'))
        silentflag = 1; % be silent
    end
end

% display banner text
RunSilent('DisplayBanner',silentflag)

% add contents to path
RunSilent('AddSubmissionContents(mfilename)',silentflag)

% Check toolboxes and versions
RunSilent('MinimumVersionChecks',silentflag)

% download required web files
RunSilent('RequiredWebFiles',silentflag)

% download required web zips
RunSilent('RequiredWebZips',silentflag)

% add contents to path (files have been downloaded)
RunSilent('AddSubmissionContents(mfilename)',silentflag)

% open an example
if ~silentflag, OpenThisFile('BrysonHo166'); end

% close this file
RunSilent('CloseThisFile(mfilename)',silentflag)

% display banner text
RunSilent('DisplayBanner',silentflag)

end
%--------------------------------------------------------------------------
function RequiredWebFiles %#ok<DEFNU>

disp('-> Obtaining required web files')

% initialize index
ind = 0;

% initialize structure
files = struct('url','','folder','');

% file 1
ind = ind + 1; % increment
files(ind).url = 'http://www1.spms.ntu.edu.sg/~lilian/bookcodes/legen/legslbdiff.m';
files(ind).folder = 'legslbdiff';

% file 2
ind = ind + 1; % increment
files(ind).url = 'http://www1.spms.ntu.edu.sg/~lilian/bookcodes/legen/lepoly.m';
files(ind).folder = 'lepoly';

% obtain full function path
full_fun_path = which(mfilename('fullpath'));
outputdir = fullfile(fileparts(full_fun_path),'include');

% download
DownloadWebFiles(files,outputdir)

disp(' ')
end
%--------------------------------------------------------------------------
function RequiredWebZips %#ok<DEFNU>

disp('-> Obtaining required web zips')

% initialize index
ind = 0;

% initialize structure
zips = struct('url','','folder','','test','');

% zip 1
ind = ind + 1; % increment
zips(ind).url = 'https://github.com/chebfun/chebfun/archive/master.zip';
zips(ind).folder = 'MFX 47023';
zips(ind).test = 'chebfun';

% zip 2
ind = ind + 1; % increment
zips(ind).url = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/8773/versions/38/download/zip/Multiprod_2009.zip';
zips(ind).folder = 'MFX 8773';
zips(ind).test = 'multiprod';

% zip 3
ind = ind + 1; % increment
zips(ind).url = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/31272/versions/10/download/zip/DataHash_20160618.zip';
zips(ind).folder = 'MFX 31272';
zips(ind).test = 'DataHash';

% zip 4
ind = ind + 1; % increment
zips(ind).url = 'http://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/40397/versions/7/download/zip/mfoldername_v2.zip';
zips(ind).folder = 'MFX 40397';
zips(ind).test = 'mfoldername';

% zip 5
ind = ind + 1; % increment
zips(ind).url = 'https://github.com/altmany/export_fig/archive/master.zip';
zips(ind).folder = 'MFX 23629';
zips(ind).test = 'export_fig';

% obtain full function path
full_fun_path = which(mfilename('fullpath'));
outputdir = fullfile(fileparts(full_fun_path),'include');

% download and unzip
DownloadWebZips(zips,outputdir)

disp(' ')
end
%--------------------------------------------------------------------------
function MinimumVersionChecks %#ok<DEFNU>

disp('-> Checking toolbox versions')

% initialize index
ind = 0;

% initialize structure
test = struct('toolbox','','version','','required','');

% test 1: MATLAB
ind = ind + 1; % increment
test(ind).toolbox = 'matlab';
test(ind).name = 'MATLAB';
test(ind).version = '0'; % any?
test(ind).required = true;

% test 2: Optimization Toolbox
ind = ind + 1; % increment
test(ind).toolbox = 'optim';
test(ind).name = 'Optimization Toolbox';
test(ind).version = '0'; % any?
test(ind).required = false;

% test 3: Parallel Computing Toolbox
ind = ind + 1; % increment
test(ind).toolbox = 'parallel';
test(ind).name = 'Parallel Computing Toolbox';
test(ind).version = '0'; % any?
test(ind).required = false;

% test 4: Symbolic Math Toolbox
ind = ind + 1; % increment
test(ind).toolbox = 'symbolic';
test(ind).name = 'Symbolic Math Toolbox';
test(ind).version = '0'; % any?
test(ind).required = false;

% test 5: Control System Toolbox
ind = ind + 1; % increment
test(ind).toolbox = 'control';
test(ind).name = 'Control System Toolbox';
test(ind).version = '0'; % any?
test(ind).required = false;

% download and unzip
VersionChecks(test)

disp(' ')
end
%--------------------------------------------------------------------------
function DisplayBanner %#ok<DEFNU>

disp('---------------------------------------------------------------')
disp('                             <strong>DTQP</strong>                              ')
disp('Link: <a href = "https://github.com/danielrherber/dt-qp-project">https://github.com/danielrherber/dt-qp-project</a>')
disp('Primary contributor: Daniel R. Herber (danielrherber on GitHub)')
disp('---------------------------------------------------------------')
end
%--------------------------------------------------------------------------
function AddSubmissionContents(name) %#ok<DEFNU>

disp('-> Adding submission contents to path')
disp(' ')

% turn off potential warning
warning('off','MATLAB:dispatcher:nameConflict')

% current file
fullfuncdir = which(name);

% current folder
submissiondir = fullfile(fileparts(fullfuncdir));

% add folders and subfolders to path
addpath(genpath(submissiondir))

% turn warning back on
warning('on','MATLAB:dispatcher:nameConflict')
end
%--------------------------------------------------------------------------
function CloseThisFile(name) %#ok<DEFNU>

disp(['-> Closing ', name])
disp(' ')

% get editor information
h = matlab.desktop.editor.getAll;

% go through all open files in the editor
for k = 1:numel(h)
    % check if this is the file
    if ~isempty(strfind(h(k).Filename,name))
        % close this file
        h(k).close
    end
end
end
%--------------------------------------------------------------------------
function OpenThisFile(name)

disp(['-> Opening ', name])

try
    % open the file
    open(name);
catch % error
    disp(['Could not open ', name])
end

disp(' ')
end
%--------------------------------------------------------------------------
function DownloadWebFiles(files,outputdir)

% store the current directory
olddir = pwd;

% create a folder for outputdir
if ~exist(outputdir, 'dir')
    mkdir(outputdir); % create the folder
else
    addpath(genpath(outputdir)); % add folders and subfolders to path
end

% change to the output directory
cd(outputdir)

% go through each file
for k = 1:length(files)

    % get data
    url = files(k).url;
    folder = files(k).folder;
    [~,nameurl,exturl] = fileparts(url);
    name = [nameurl,exturl];

    % first check if the test file is in the path
    if exist(name,'file') == 0

        try
            % download file
            outfilename = websave(name,url);

            % create a folder utilizing name as the foldername name
            if ~exist(fullfile(outputdir,folder), 'dir')
                mkdir(fullfile(outputdir,folder));
            end

            % move the file
            movefile(outfilename,fullfile(outputdir,folder))

            % output to the command window
            disp(['Downloaded ',folder,'/',name])

        catch % failed to download
            % output to the command window
            disp(['Failed to download ',folder,'/',name])

            % remove the html file
            delete([name,'.html'])
        end

    else
        % output to the command window
        disp(['Already available ',name])
    end
end

% change back to the original directory
cd(olddir)
end
%--------------------------------------------------------------------------
function DownloadWebZips(zips,outputdir)

% store the current directory
olddir = pwd;

% create a folder for outputdir
if ~exist(outputdir, 'dir')
    mkdir(outputdir); % create the folder
else
    addpath(genpath(outputdir)); % add folders and subfolders to path
end

% change to the output directory
cd(outputdir)

% go through each zip
for k = 1:length(zips)

    % get data
    url = zips(k).url;
    folder = zips(k).folder;
    test = zips(k).test;

    % first check if the test file is in the path
    if exist(test,'file') == 0

        try
            % download zip file
            zipname = websave(folder,url);

            % save location
            outputdirname = fullfile(outputdir,folder);

            % create a folder utilizing name as the foldername name
            if ~exist(outputdirname, 'dir')
                mkdir(outputdirname);
            end

            % unzip the zip
            unzip(zipname,outputdirname);

            % delete the zip file
            delete([folder,'.zip'])

            % output to the command window
            disp(['Downloaded and unzipped ',folder])

        catch % failed to download
            % output to the command window
            disp(['Failed to download ',folder])

            % remove the html file
            delete([folder,'.html'])
        end

    else
        % output to the command window
        disp(['Already available ',folder])
    end
end

% change back to the original directory
cd(olddir)
end
%--------------------------------------------------------------------------
function RunSilent(str,silentflag)

% if silent, capture the output
if silentflag
    O = evalc(str); %#ok<NASGU>
else
    eval(str);
end
end
%--------------------------------------------------------------------------
function VersionChecks(test)

% initialize counter
counter = 0;

% go through each file
for k = 1:length(test)
    try
        if verLessThan(test(k).toolbox,test(k).version) % failed
            if test(k).required % required
                str = ['Failed (REQUIRED): ',test(k).name];
            else % recommended
                str = ['Failed (optional): ',test(k).name];
            end

        else % passed
            str = ['Passed: ',test(k).name];
            counter = counter + 1;
        end

    catch % failed to check the toolbox
        str = ['Failed to check toolbox: ', test(k).name];

    end
    if ~strcmpi(test(k).version,'0')
        str = [str,' -v', test(k).version]; %#ok<AGROW>
    end

    % display to command window
    disp(str)
end

% check if all tests were passed
if counter == length(test) % successful
    disp('All toolbox and version checks passed')
else % failure
    warning('Not all toolbox and version checks were successful')
end
end