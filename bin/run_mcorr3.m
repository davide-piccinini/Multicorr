function [] = run_mcorr(configfile)

%% 25/11/2023 aggiunto il reset del calcolo parallelo ogni 25 giorni

delete(gcp('nocreate'))

fprintf('Process started at %s\n',datestr(now,'yyyy-mm-ddTHH:MM:SS'))
YTIME=datetime('now');
% clear all



%% Read Configuration File
if nargin < 1
    CFG.FULLDAY = {};
    CFG.ROOT    = [];
    CFG.DAYDIR  = './';
    CFG.FILT    = [2 15];
    CFG.CCTHR   = 0.55;
    CFG.NP      = 3;
    CFG.NS      = 3;
    CFG.PLWIN   = 0.7; % P wave template length
    CFG.SLWIN   = 1.0; % S wave template length
    CFG.PREWIN  = 0.2; % Add this winlen before the templates
    CFG.CORRECT = 1;     %% FIXED (0) or ADAPTIVE(1) P template length
    CFG.TEMPLATESDIR = 'TEMPLATES/';
    CFG.WORKDIR      = 'WORK/';
    CFG.MATDIR  = 'OUT/';
    CFG.NLLOUT  = 'NLL_OBS/';
    CFG.DT      = 5; % minimum distance (s) between two matches
    CFG.FIG     = 1; % not used
    CFG.GARBAGE = 1; % clear WORK directory
    CFG.NWORK   = 16;% number of cores used
    CFG.SEARCH  = -1; % scan only for templates occurred N days after and before; set -1 for no selection
    CFG.NMIN    = 10 ; % minimum number of events in case of no events in CFG.SEARCH selection
    CFG.REFRESH = 0; % restart parpool each day to save memory (0/1)
    CFG.ADDPICK = 0; % add picks from sac files (dir SAC in TEMPLATES)
    CFG.JSON    = 0; % read json formatted file as template
else
    CFG=readcfg(configfile);  %% Now use the readcfg.m function (added 5 oct 2022)
end

NF=fieldnames(CFG);
fprintf('Running using following parameters:\n')
O=fopen('Running_parameters.txt','a');
fprintf(O,'%s\n',datestr(now));
for k=1:length(NF);
    VAR=getfield(CFG,char(NF(k)));
    fprintf('%s = %s\n',char(NF(k)),string(VAR))
    fprintf(O,'%s = %s\n',char(NF(k)),string(VAR));
end
fclose(O);
%% Create OUT directory
[SUCCESS,MESSAGE,MESSAGEID] = mkdir('OUT');
if SUCCESS ~=1
MESSAGE
end
%% Create template list
if CFG.JSON==0
D=dir([CFG.TEMPLATESDIR '*.xml']);
fout=fopen('lista','w');
for j=1:length(D)
    CD=char(D(j).name);
    fprintf(fout,'%s\n',CD(1:end-4));
end
fclose(fout);
else
D=dir([CFG.TEMPLATESDIR '*.json']);
fout=fopen('lista','w');
for j=1:length(D)
    CD=char(D(j).name);
    fprintf(fout,'%s\n',CD(1:end-5));
end
fclose(fout);
end    
%% Read in template list
A=textread('lista','%s');
idA=[]; %% put here templates id to be removed
A(idA)=[];
%% or 
%idA=[]; % put here DDIR id to 
%A=A(idA);
%%
fprintf('Number of Templates: %3d  -- %s -> %s\n',length(A),char(A(1)),char(A(end)))

for k=1:length(A);
    MA(k)=datenum(A{k},'yyyy-mm-ddTHH:MM:SS.FFF');
end
%keyboard

%% Read daylong directories
DDIR = dir([CFG.DAYDIR '20*']);
%idD=[]; % put here DDIR id to be removed
%DDIR(idD)=[];
% otherwise put here directory to process
idD=[1:2]; % put here DDIR id to process
DDIR=DDIR(idD);
fprintf('Number of Days: %3d -- %s -> %s\n',length(DDIR),DDIR(1).name,DDIR(end).name)

%%
%% check if parpool is running (check this! it seems that something goes wrong when invoked inside a function)
% if isempty(gcp('nocreate'))
%     mypool=parpool
% end

%% build a superstructure with all the days
FULLSTRUCT={};
try
    fprintf('Wake up Parallel pool for the first time\n')
    p=parpool(CFG.NWORK);
catch
    while exist('p') ==0 % p.Connected==0;
        fprintf('It didnt wake up, try again\n');
        p=parpool(CFG.NWORK);
    end
end
for i=1:numel(DDIR)
    FULLSTRUCT{i}=mcorr_io_readDAY(['./' DDIR(i).name] ,CFG.FILT); 
end
if CFG.REFRESH==1
    fprintf('---Closing day reading parpool---\n')
    delete(p); % close read Day parpool object 
end
%% INITIALIZE PARPOOL IL REFRESH==0
if CFG.REFRESH==0;
    try
        p=parpool(CFG.NWORK);
    catch
        while exist('p')==0 % p.Connected==0;
            fprintf('Trying to wake up parallel pool\n');
            p=parpool(CFG.NWORK);
        end
    end
end

%% ------------------- BEGIN MAIN CYCLE ----------------
for j=1:numel(FULLSTRUCT)
    p = gcp('nocreate');
    if exist('p')==0     %% se è 0 vuol dire che non è stato creato prima -> REFRESH=1     
        try
            p=parpool(CFG.NWORK);
        catch
            while exist('p')==0 % p.Connected==0;
                fprintf('Trying to wake up parallel pool\n');
                p=parpool(CFG.NWORK);
            end
        end
    end
  % keyboard  
    DAYDIR =FULLSTRUCT{j}{1}.DAYDIR; DAYDIR=DAYDIR(3:end);
    MDAYDIR=datenum(DAYDIR,'yyyymmdd');
    
    %% TEMPLATES SELECTION ----------------------------------
    if CFG.SEARCH > 0 % TEMPLATES SELECTION TEMPORAL OR RANDOM
        idx=find(MA > MDAYDIR-CFG.SEARCH & MA <= MDAYDIR+CFG.SEARCH+1); 
        fprintf('Selecting a %2.0f days pool: %d\n',2*CFG.SEARCH,numel(idx))
        
        if numel(idx)>= CFG.NMIN
            TA=A(idx);
            fprintf('Number of templates used: %d\n\n',numel(idx))
        end
        if numel(idx)==0
           fprintf('Selecting a %2.0f full random pool\n',CFG.NMIN)
           idx=randperm(numel(A),CFG.NMIN);
           TA=A(idx);
           fprintf('Number of templates used: %d\n\n',numel(idx))           
        end
        if numel(idx) > 0 & numel(idx) < CFG.NMIN
           fprintf('Adding a %2.0f partial random pool\n',CFG.NMIN)
           idx=[idx, randperm(numel(A),CFG.NMIN)];
           TA=A(idx);
           fprintf('Number of templates used: %d\n\n',numel(idx))
        end
    elseif CFG.SEARCH == 0 %% TEMPLATES ARE A RANDOM subset OF THE TEMPLATES OCCURRED DURING THE DAY
        index=find(MA > MDAYDIR-CFG.SEARCH & MA <= MDAYDIR+CFG.SEARCH+1);
        if numel(index) >= CFG.NMIN;
            idx=randperm(numel(index),CFG.NMIN);
        else            
            iADD=0;
            while numel(index) < CFG.NMIN
                iADD=iADD+1;
                index=find(MA > MDAYDIR-CFG.SEARCH -iADD & MA <= MDAYDIR+CFG.SEARCH+1+iADD);
            end            
            idx=randperm(numel(index),CFG.NMIN);
        end
        
        TA=A(index(idx));
        fprintf('Number of templates used: %d\n\n',numel(idx))
    else  %% NO SELECTION --> USE FULL SET OF TEMPLATES
        TA=A;
    end
    %% ----------------------------------------------------------------
    
    FULLDAY=FULLSTRUCT{j};
    DAYDIR =FULLDAY{1}.DAYDIR; DAYDIR=DAYDIR(3:end);
%    keyboard
    for k=1:length(TA)
        LA=char(TA(k));
        root=LA;
        fprintf('%s -- Template n. %3.0f / %3.0f - Day n. %3.0f / %3.0f\n',DAYDIR, k,length(TA),j,numel(FULLSTRUCT))        
        CFG.ROOT=root;
        CFG.DAYDIR=DAYDIR;
        CFG.FULLDAY=FULLDAY;
        mcorr ('FULLDAY', CFG.FULLDAY,...  
            'ROOT', CFG.ROOT,...
            'DAYDIR', CFG.DAYDIR,...
            'TEMPLATESDIR', CFG.TEMPLATESDIR, ...
            'WORKDIR', CFG.WORKDIR, ...
            'PLWIN', CFG.PLWIN, ...
            'SLWIN', CFG.SLWIN, ...
            'PREWIN', CFG.PREWIN, ...
            'FILT', CFG.FILT, ...
            'CORRECT', CFG.CORRECT, ...
            'NP', CFG.NP, ...
            'NS', CFG.NS, ...
            'CCTHR', CFG.CCTHR, ...
            'GARBAGE', CFG.GARBAGE,...
             'ADDPICK',CFG.ADDPICK,...
             'JSON',CFG.JSON);
    end
    % pack %% inserito il 25 novembre vediamo se migliora la gestione della memoria
    if CFG.REFRESH==1
        poolobj = gcp('nocreate');
        delete(poolobj)
    end
    if ismember(j,10:10:200); % every X days restart parallelpool
        poolobj = gcp('nocreate');
        delete(poolobj);
        
        p = gcp('nocreate');
        if numel(p)==0     %% se è 0 vuol dire che non è stato creato prima -> REFRESH=1
            try
                fprintf('Trying to run parallel pool\n');
                p=parpool(CFG.NWORK);
            catch
                while numel(p)==0 % p.Connected==0;
                    fprintf('No way..!! Trying to wake up parallel pool\n');
                    p=parpool(CFG.NWORK);
                end
            end
        end
        
    end
end

delete (gcp)

ETIME=datetime('now');
DURATION=ETIME-YTIME;
fprintf('Total Time Used:  %s\n',DURATION)
fprintf('Process ended at %s\n',datestr(now,'yyyy-mm-ddTHH:MM:SS'))




