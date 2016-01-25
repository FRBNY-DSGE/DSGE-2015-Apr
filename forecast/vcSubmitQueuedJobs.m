function JobOut = vcSubmitQueuedJobs(nJobs,JobFcnHandle,nJobArgOut,JobOptions,...
    ListPathDependencies,nMaxWorkers,varargin)

% vcSubmitQueuedJobs
%
% Submit parallel jobs that might exceed the maximum number of workers
% available.
%
% Usage:
%   vcSubmitQueuedJobs(nJobs,JobFcnHandle,nJobArgOut,JobOptions,ListPathDependencies)
%   vcSubmitQueuedJobs(...,OptionName,OptionValue,...)
%
% Mandatory Inputs:
%   
%   nJobs
%   Number of Jobs to run
%
%   JobFcnHandle
%   Handle to the function called by each worker.
%
%   nJobArgOut
%   Number of output arguments in each job.
%   
%   JobOptions
%   Cell array with options required by the function called by the worker.
%   It needs to contain one element per job, unless no inputs are
%   necessary, in which case, {} needs to be used.
%
%   ListPathDependencies
%   Cell array with the list of path dependencies needed in function called
%   by workers.
%
%   nMaxWorkers
%   Maximum number of simultaneous workers. 
%
% Optional inputs:
%
%   nMaxTasks
%   Maximum number of tasks per job. Default: 1.
%
%   ShowOutput
%   If set to 1 all output of each worker is shown. 
%   Default: 1 if LogFileName is empty, 0 otherwise. 
%   Note: even though the default is set depending on LogFileName, this
%   option can be set independendently of LogFileName.
%
%   LogFileName
%   If specified, then all output is saved in a ascii file.
%
%   OutputPreamble
%   If specified this should be a string to introduce in a fprintf command
%   before showing the output of each worker. It must receive as input the
%   worker number.
%
%   SaveTmpFileName
%   If specified, it will save the environment in a mat file with the name 
%   provided. This file is deleted upon completion and data gathering from
%   all jobs.
%
%   SaveTmpVarAll
%   If set to 1 it will save the whole environment, otherwise it will save 
%   only selected information.
%   Default: 0;
%   NOTE: Only relevant if SaveTmpFileName is specified as well.
%
%   SaveTmpDelete
%   If set to one the temporary file is deleted on conclusion. If set to 
%   zero the temporary file is left in the folder.
%   Default: 1
%
% .........................................................................
%
% Created: February 11, 2009
% Updated: August 24, 2010
% by Vasco Curdia

%% ------------------------------------------------------------------------

%% default options
nMaxTasks = 1;
SaveTmpDelete = 1;
SaveTmpVarAll = 0;
LogFileName = '';
PauseTime = 10;

%% Update options
if ~isempty(varargin)
    nOptions = length(varargin);
    if mod(nOptions,2), error('Incorrect number of optional arguments.'), end
    for jO=1:nOptions/2
        eval(sprintf('%s = varargin{%.0f};',varargin{(jO-1)*2+1},jO*2))
    end
end

%% Checks
if isempty(JobOptions),for j=1:nJobs,JobOptions{j} = {};end,end
if ~exist('ShowOutput','var'),ShowOutput = isempty(LogFileName);end
if ~exist('SaveTmpVarList','var')
    SaveTmpVarList = {'nJobs','JobFcnHandle','nJobArgOut','JobOptions',...
        'ListPathDependencies','nMaxWorkers','varargin','nCompleted',...
        'nSubmitted','j','jJob','jobsRunning','JobOut','JobTasks','JobIdx'};
end

%% Submit jobs
nCompleted = 0;
nSubmitted = 0;
j = 0;
jJob = 0;
jobsRunning = false(nJobs,1);
JobOut = cell(nJobs,nJobArgOut);
JobTasks=[];
JobIdx=[];
%JobToc = toc;
if exist('SaveTmpFileName','var')
    if SaveTmpVarAll
        save(SaveTmpFileName)
    else
        save(SaveTmpFileName,SaveTmpVarList{:})
    end
end
while nCompleted<nJobs
    while (nSubmitted<nMaxWorkers) && (j<nJobs)
        jJob = jJob+1;
        job{jJob} = createSgeJob(0);
        set(job{jJob},'PathDependencies',ListPathDependencies)
        JobIdx{jJob} = [];
        for jTask=1:min([nMaxTasks,nJobs-j,nMaxWorkers-nSubmitted])
            j = j+1;
            JobIdx{jJob}(end+1) = j;
            nSubmitted = nSubmitted+1;
            TaskID{jJob,jTask} = createTask(job{jJob},JobFcnHandle,nJobArgOut,JobOptions{j});
            if ShowOutput || ~isempty(LogFileName)
                set(TaskID{jJob,jTask}, 'CaptureCommandWindowOutput', true);
            end
        end
        JobTasks(jJob) = jTask;
        submit(job{jJob})
        jobsRunning(jJob) = true;
    end
    pause(PauseTime)
    % save tmp file if desired
    if exist('SaveTmpFileName','var')
        if SaveTmpVarAll
            save(SaveTmpFileName)
        else
            save(SaveTmpFileName,SaveTmpVarList{:})
        end
    end
    listJobsRunning = find(jobsRunning);
    for jj=1:length(listJobsRunning)
        jobj = listJobsRunning(jj);
%        fprintf('%s: job %.0f, %s\n',vctoc(JobToc),jobj,get(job{jobj},'State'))
        % save log, if desired
        if ~isempty(LogFileName)
            for jTask=1:JobTasks(jobj)
                jobLogName = sprintf('%s%03.0f.log',LogFileName,JobIdx{jobj}(jTask));
                fid = fopen(jobLogName,'wt');
                if exist('OutputPreamble','var')
                    fprintf(fid,OutputPreamble,JobIdx{jobj}(jTask));
                end
                fprintf(fid,strrep(get(TaskID{jobj,jTask},'CommandWindowOutput'),'%','%%'));
                if ~isempty(get(TaskID{jobj,jTask},'errormessage'))
                    fprintf(fid,strrep(get(TaskID{jobj,jTask},'erroridentifier'),'%','%%'));fprintf('\n');
                    fprintf(fid,strrep(get(TaskID{jobj,jTask},'errormessage'),'%','%%'));fprintf('\n');
                end
                fclose(fid);
            end
        end
        if strcmp(get(job{jobj},'State'),'finished')
            jobsRunning(jobj) = false;
            nSubmitted = nSubmitted-JobTasks(jobj);
            nCompleted = nCompleted+JobTasks(jobj);
            % collect output
            if ~isempty(getAllOutputArguments(job{jobj}))
                JobOut(JobIdx{jobj},:) = getAllOutputArguments(job{jobj});
            else
                fprintf('Note: Output for Job %i is empty\n',jobj);
            end
            % save log, if desired
            if ~isempty(LogFileName)
                for jTask=1:JobTasks(jobj)
                    jobLogName = sprintf('%s%03.0f.log',LogFileName,JobIdx{jobj}(jTask));
                    fid = fopen(jobLogName,'wt');
                    if exist('OutputPreamble','var')
                        fprintf(fid,OutputPreamble,JobIdx{jobj}(jTask));
                    end
                    fprintf(fid,strrep(get(TaskID{jobj,jTask},'CommandWindowOutput'),'%','%%'));
                    if ~isempty(get(TaskID{jobj,jTask},'errormessage'))
                        fprintf(fid,strrep(get(TaskID{jobj,jTask},'erroridentifier'),'%','%%'));fprintf('\n');
                        fprintf(fid,strrep(get(TaskID{jobj,jTask},'errormessage'),'%','%%'));fprintf('\n');
                    end
                    fclose(fid);
                end
            end
%             % save tmp file if desired
%             if exist('SaveTmpFileName','var')
%                 if SaveTmpVarAll
%                     save(SaveTmpFileName)
%                 else
%                     save(SaveTmpFileName,SaveTmpVarList{:})
%                 end
%             end
        end
    end
end
% save tmp file if desired
if exist('SaveTmpFileName','var')
    if SaveTmpVarAll
        save(SaveTmpFileName)
    else
        save(SaveTmpFileName,SaveTmpVarList{:})
    end
end
nJob = jJob;
j=0;
for jJob=1:nJob
    for jTask=1:JobTasks(jJob)
        j=j+1;
        % show onscreen output if desired
        if ShowOutput
            if exist('OutputPreamble','var')
                fprintf(OutputPreamble,j);
            end
            fprintf(strrep(get(TaskID{jJob,jTask},'CommandWindowOutput'),'%','%%'));
        end
        if ~isempty(get(TaskID{jJob,jTask},'errormessage'))
            fprintf(strrep(get(TaskID{jJob,jTask},'erroridentifier'),'%','%%'));fprintf('\n');
            fprintf(strrep(get(TaskID{jJob,jTask},'errormessage'),'%','%%'));fprintf('\n');
        end
    end
    % kill job
    destroy(job{jJob})
end
% delete tmp file
if exist('SaveTmpFileName','var') && SaveTmpDelete
    if strcmp(SaveTmpFileName(max(end-3,1):end),'.mat')
        delete(SaveTmpFileName)
    else
        delete([SaveTmpFileName,'.mat'])
    end
end


