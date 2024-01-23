function mcorr(varargin)

% Usage: function mcorr(rootxml,roottemp,DAYDIR)

% Multichannel xcorr
% dipendenze: mcorr_io_ReadMSEED
%             mcorr_core_findwavmseed.m
%
% Legge QuAKEML catalog e i relativi files mseed
% Estrae i template e fa la crosscorrelazione
% Scrive i risultati di ogni canale/fase nella dir WORK
% Al termine crea una struttura matlab nella dir OUT
% GLI EVENTI DEVONO AVERE PICK P e S
%
% per preparare i templates usare il codice python alla fine del codice
%
% -----------------------------------------
% usage:    mcorr(rootxml,roottemp,DAYDIR)
% -----------------------------------------
%
% Davide Piccinini Dec 2020 mailto:davide.piccinini@ingv.it
% mcorr ver. 1.1

% 04/02/2021
% + Velocizzata l'associazione degli eventi
% + cambiato il modo di lettura dei dati mseed (da mcorr_io_ReadMSEEDFast a mcorr_io_rdmseed)
% + incluse tutte le funzioni associate (findwav e mcorr_io_rdmseed)
% 17/02/2021
% + aggiunto in mcorr_core_findwavmseed il controllo sulla maxCC dei doppioni
% + aggiunge il valore di PREWIN ai tempi degli slaves
% 18/02/2021
% + aggiunto un controllo sulla PREF e SREF
% 17/12/2021
% + upgraded mcorr_io_decodeQML to correctly get only real picks
% 06/08/2022
% + modificato il ciclo if exist('IS') ==1 precedentemente era if IS > 0
% 25/06/2022
% + modified mcorr_io_decodeQML
% 7/10/2022 
% + changed VAR name wrong CCTHR in PLAG & SLAG at the end of mcorr
% 12/10/2022
% + deeply modified mcorr_core_findwaw call 
% + Added waveforms amplitude (in counts) of the matches and location of the template in the .mat output;
% + USE IT to estimate magnitude (to be done!!!)  
% + added also P and S maximum of the template
% 11/02/2023
% Added waveforms for both template and match in the SLAVE output 
% Added external Picks input (read SAC files with P and S time)
% Added WEIGHT of the template (to be used with selectcorr)

%% INPUTS

Def=struct( ...
    'FULLDAY',{},...
    'ROOT',[],...
    'DAYDIR',[],...
    'TEMPLATESDIR', [], ...
    'WORKDIR',[], ... 
    'PLWIN',0.5, ...          % sec of P
    'SLWIN',0.8, ...          % sec of S
    'PREWIN', 0.0, ...        % add this before P and S
    'CORRECT', 1, ...         % time correction for PLWIN sets PLWIN eq to S-P
    'FILT', [2 15], ...       % BandPass filter extrema
    'NP', 4, ...              % Numero minimo di fasi P
    'NS', 3, ...              % numero minimo di fasi S    
    'CCTHR', 0.65, ...        % soglia di CC
    'GARBAGE', 1,...          % clear WORK directory each template
    'ADDPICK',0, ...          % get picks from outside
    'JSON',0 ...              % use JSON instead !QML
    );
Args=mcorr_scaffold_parseargs(Def,varargin);

%
if isempty(Args.TEMPLATESDIR) | isempty(Args.TEMPLATESDIR)
    error(['### mcorr.m: ERROR !! Please Specify a templates-dirpath ' ...
           '(TEMPLATESDIR) and a working-dirpath (WORKDIR)!!'])
end

%% SETUP WORKING DIR
if exist(Args.WORKDIR,'dir')==0
    mkdir(Args.WORKDIR)
end

if exist(Args.TEMPLATESDIR,'dir')==0
    error('### mcorr.m: ERROR !! TEMPLATESDIR doesn''t exists!!')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LEGGE L'EVENTO
if Args.JSON==0
    fprintf('Getting template picks from %s\n',Args.ROOT)
    filexml=[Args.TEMPLATESDIR '/' Args.ROOT '.xml'];
    X=mcorr_io_decodeQML(filexml,0);
end
if Args.JSON==1
    fprintf('Getting template picks from %s\n',Args.ROOT)
    filexml=[Args.TEMPLATESDIR '/' Args.ROOT '.json'];
    X=mcorr_io_decodeJSON(filexml);
end
    
 %keyboard

%% esce se il n di fasi P ed S non è sufficiente
namex=fieldnames(X);
np=0; ns=0;
for k=1:length(namex)-1
    XU=getfield(X,char(namex(k)));
    if XU.P > 0; np=np+1; end
    if XU.S > 0; ns=ns+1; end
end
if np < Args.NP || ns < Args.NS ; fprintf('Template has insufficient number of P and S [%d %d / %d %d]\n',np,ns,Args.NP,Args.NS);return; end

%% LEGGE IL TEMPLATE MSEED
D=dir([Args.TEMPLATESDIR '/' Args.ROOT '.*.ms']);
if isempty(D)==1
    fprintf('No miniseed data for template\n');
    return
end
nn=0;
for ii=1:length(D)
    nn=nn+1;
    S(nn)=mcorr_io_ReadMSEED([D(ii).folder '/' D(ii).name ]);
end

for j=1:length(S)
    S(j).station=strrep(S(j).station,' ',''); %% accounting for 3char station code
    trace=S(j).data;
    trace=trace-mean(trace);
    trace=mcorr_signal_taper(trace,50);
    S(j).datafilt = mcorr_signal_filtrax(trace,Args.FILT(1),Args.FILT(2),S(j).sampleRate);
end

n=0;
FN=fieldnames(X);
for k=1:length(S)
    NM=S(k).station;
%    NM=strrep(NM,' ',''); %% accounting for 3char station code
    NM=deblank(NM); %% accounting for 3char station code
    for jj=1:length(FN)-1
        if strcmp(char(FN(jj)),NM)==1
            PT=X.(NM).P;
            Punc=X.(NM).Punc;
            ST=X.(NM).S;
            Sunc=X.(NM).Sunc;
            if Args.CORRECT==1 && isempty (PT) ==0 && isempty (ST) ==0
                %% check P win and resize to SP time
                SPT=86400*(ST-PT);
                if Args.PLWIN > SPT
                    Args.PLWIN=SPT;
                    fprintf('Resizing Pwin length @%s to %4.2f s\n',NM,Args.PLWIN)
                end
            end
            if isempty (PT) == 0
                n=n+1;
                [~, ps]=min(abs(S(k).matlabTimeVector-PT));
                S(k).PTIME=PT;
                S(k).Punc=Punc;
                try
                    S(k).ptemplate=(S(k).datafilt(round(ps-(Args.PREWIN*S(k).sampleRate)):round(ps+(Args.PLWIN*S(k).sampleRate))));
                    if  ST > 0
                        [~, ss]=min(abs(S(k).matlabTimeVector-ST));
                        S(k).STIME=ST;
                        S(k).Sunc=Sunc;
                        S(k).stemplate=(S(k).datafilt(ss-(Args.PREWIN*S(k).sampleRate):ss+(Args.SLWIN*S(k).sampleRate)));
                    else
                        S(k).STIME=[];
                        S(k).Sunc =[];
                        S(k).stemplate=[];
                    end
                    NS(n)=S(k);
                catch
                    fprintf('Something goes wrong with the pick [Maybe S outside the template?] @%s\n',NM)
                    n=n-1;
                end
            end
        end
    end
end

% keyboard

%% HERE CALL THE FUNCTION TO ADD MANUAL PICKS 
if Args.ADDPICK==1
    LNS=length(NS);
    % SCAN THE TEMPLATES to check if there are some sac files
    DSAC= dir([Args.TEMPLATESDIR '/SAC/' Args.ROOT(1:end-4) '*']);  
    fprintf('---> ADDING %3d sac files for this event\n',length(DSAC))
    for ksac=1:length(DSAC)  
        PT=[];ST=[];
        PICK=mcorr_io_get_picks([Args.TEMPLATESDIR '/SAC/' DSAC(ksac).name]);
        NS(LNS+ksac).network=PICK.net;
        NS(LNS+ksac).station=PICK.sta;
        NS(LNS+ksac).sampleRate=round(1/PICK.DT);
        NS(LNS+ksac).matlabTimeVector=datenum(PICK.time);
        NS(LNS+ksac).data=PICK.data;

        trace=PICK.data;
        trace=trace-mean(trace);
        trace=mcorr_signal_taper(trace,50);
        NS(LNS+ksac).datafilt = mcorr_signal_filtrax(trace,Args.FILT(1),Args.FILT(2),NS(ksac).sampleRate);
        NS(LNS+ksac).PTIME=datenum(PICK.P);
        NS(LNS+ksac).Punc=datenum(PICK.Punc);
        NS(LNS+ksac).STIME=datenum(PICK.S);
        NS(LNS+ksac).Sunc=datenum(PICK.Sunc);
        NS(LNS+ksac).channel=PICK.CHN;
        PT=NS(LNS+ksac).PTIME;ST=NS(LNS+ksac).STIME;
        % getting template wvf
        if Args.CORRECT==1 && isempty (PT) ==0 && isempty (ST) ==0
            %% check P win and resize to SP time
            SPT=86400*(NS(LNS+ksac).STIME-NS(LNS+ksac).PTIME);
            if Args.PLWIN > SPT
                Args.PLWIN=SPT;
                fprintf('Resizing Pwin length @%s to %4.2f s\n',NS(LNS+ksac).station,Args.PLWIN)
            end
        end

        if isempty (PT) == 0
            n=n+1;
            [~, ps]=min(abs(NS(LNS+ksac).matlabTimeVector-PT));
            try
                NS(LNS+ksac).ptemplate=(NS(LNS+ksac).datafilt(round(ps-(Args.PREWIN*NS(LNS+ksac).sampleRate)):round(ps+(Args.PLWIN*NS(LNS+ksac).sampleRate))));
                if  ST > 0
                    [~, ss]=min(abs(NS(LNS+ksac).matlabTimeVector-ST));
                    NS(LNS+ksac).stemplate=(NS(LNS+ksac).datafilt(ss-(Args.PREWIN*NS(LNS+ksac).sampleRate):ss+(Args.SLWIN*NS(LNS+ksac).sampleRate)));
                end
            catch
                fprintf('Something goes wrong with the SAC pick [Maybe S outside the template?] @%s\n',NS(LNS+ksac).station)
                n=n-1;
            end
        end
    end
end
%% END OF SAC PICK 



%% CHECK AGAIN P AND S NUMBER
np=0;
ns=0;
for k=1:length(NS)
    if NS(k).PTIME > 0
        np=np+1;
    end
    if NS(k).STIME > 0
        ns=ns+1;
    end
end
if np < Args.NP || ns < Args.NS; fprintf('Insufficient phases --- STOP ---\n'); return; end


for k=1:length(NS)   %% ADDING ORIGIN TIME
    NS(k).OTime=char(X.HEADER(1));
    NS(k).MOTime=cell2mat(X.HEADER(2));
    NS(k).LAT=cell2mat(X.HEADER(3));
    NS(k).LON=cell2mat(X.HEADER(4));
    NS(k).DEP=cell2mat(X.HEADER(5));
    NS(k).MAG=cell2mat(X.HEADER(6));
end
OTime=NS(1).OTime;

fprintf('Number of stations : %3.0f (%3.0f P %3.0f S)\n',length(NS),np,ns)

T0=clock;
fprintf('Correlating template vs continuous...\n')
%keyboard
FULLDAY=Args.FULLDAY;
CCTHR=Args.CCTHR;

Args.FULLDAY=[];
varargin{1}=[];
varargin{2}=[];

%% mod 12/10/2022 
n=1;
FULLDAYIDx={};
NSID=[];
for k=1:numel(NS);
    for j=1:numel(FULLDAY)
        if strcmp(deblank(NS(k).station),deblank(FULLDAY{j}.station))==1
            FULLDAYIDx{n}=FULLDAY{j}; %% queste sono solo le stazioni presenti nel template
            NSID(n)=k;
            n=n+1;            
        end
    end
end
NS=NS(NSID); %% rimuove stazioni che sono nel template ma non nel continuo
%% ALL'USCITA del ciclo NS e FULLDAYIDx hanno stesso ordine e stesso numero di stazioni
%for ix=1:numel(NS)
%parfor ix=1:numel(NS)
%    mcorr_core_findwavmseed(FULLDAYIDx{ix},NS(ix),CCTHR,'P');
%    if isempty(getfield(NS(ix),'STIME'))==0
%        mcorr_core_findwavmseed(FULLDAYIDx{ix},NS(ix),CCTHR,'S');
%    end
%end

% keyboard
np=0;
ns=0;
for ix=1:numel(NS)
	np=np+1;
    FP(np)=parfeval(@mcorr_core_findwavmseed,0,FULLDAYIDx{ix},NS(ix),CCTHR,'P');
    if isempty(getfield(NS(ix),'STIME'))==0
	    ns=ns+1;
         FS(ns)=parfeval(@mcorr_core_findwavmseed,0,FULLDAYIDx{ix},NS(ix),CCTHR,'S');
    end
end
try
wait(FP);
FP;
wait(FS);
FS;
catch
	fprintf('Parfeval wait fail...<n')
end
clear FP FS 

 parfevalOnAll(@clearvars,0);

 % clear FULLDAYIDx NSID S 

%% DEBUG 
%KK=whos;
%save(['DEBUG/' OTime '.' Args.DAYDIR '.memory.mat'], 'KK');

% keyboard
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ET=etime(clock,T0);
fprintf('Total time used for correlation: %8.2f sec\n',ET)
%keyboard
%% %%%%%%%%%%%%%%%%%%%%
%% ASSOCIA LE P E LE S
%% ====================
T0=clock; 
try
T = struct2table(NS);
SS=sortrows(T,'PTIME');
NS=table2struct(SS);
clear SS
catch
    fprintf('NS is one STOP ********** \n');
    return
end

for k=1:length(NS);  % uncomment for non parallel
    STA=NS(k).station;
    PTEMP(k)=NS(k).PTIME;
    try        
        Dp=dir([Args.WORKDIR OTime '.' STA '*.P.' Args.DAYDIR '.out']);
        if Dp.bytes ~= 0
        % fprintf('Reading %s\n',[ Args.WORKDIR OTime '.' STA '*.P.' Args.DAYDIR '.out'])
            [PT, CORR(k).AP, CORR(k).CCP]=textread([Args.WORKDIR Dp(1).name],'%s %f %f');
            CORR(k).PT=datenum(PT,'yyyy-mm-ddTHH:MM:SS.FFF');
        else
            fprintf('No P time for Station %s\n',STA);
            CORR(k).AP=NaN;
            CORR(k).CCP=NaN;
            CORR(k).PT=NaN;
            PTEMP(k)=NaN;
        end

    catch ME
        fprintf('No file P time for Station %s\n',STA);
        CORR(k).AP=NaN;
        CORR(k).CCP=NaN;
        CORR(k).PT=NaN;
        PTEMP(k)=NaN;
    end
    if k==1 & isempty(CORR(k).AP)==1; fprintf('No P-events found for first station\n');if Args.GARBAGE==1; fprintf('purging working directory...\n');delete('WORK/*');end ;clear all; return; end  %% SE NON TROVA P ESCE


    try
        STEMP(k)=NS(k).STIME;
        Ds=dir([Args.WORKDIR OTime '.' STA '*.S.' Args.DAYDIR '.out']);
        % fprintf('Reading %s\n',['WORK/' OTime '.' STA '*.S.' DAYDIR '.out'])
        if Ds.bytes ~=0
            [ST, CORR(k).AS, CORR(k).CCS]=textread([Args.WORKDIR Ds(1).name],'%s %f %f');
            CORR(k).ST=datenum(ST,'yyyy-mm-ddTHH:MM:SS.FFF');
        else
            fprintf('No S time for Station %s\n',STA);
            STEMP(k)=NaN;
            CORR(k).AS=NaN;
            CORR(k).CCS=NaN;
            CORR(k).ST=NaN;
        end

    catch ME
        fprintf('No file S time for Station %s\n',STA);
        STEMP(k)=NaN;
        CORR(k).AS=NaN;
        CORR(k).CCS=NaN;
        CORR(k).ST=NaN;
    end
end

%% STABILISCE LA STAZIONE DI RIFERIMENTO PER IL CALCOLO DELLE TT
SEC=1/86400;
PREF=1; %STAZIONE DI RIFERIMENTO PER LE P E' LA PRIMA STAZIONE
if isnan(CORR(PREF).CCP)==1
    for k=1:length(CORR)
        n(k)=length(CORR(k).CCP);
    end
    is=find(n>1);
    PREF=min(is);
end

%% SREF POTREBBE NON ESSERE LA 1° STAZ PERCHE' POTREBBE NON AVERE LA P
for k=1:length(NS)
    if isempty(NS(k).STIME)
        SREF(k)=0;
    else
        SREF(k)=NS(k).STIME;
    end
end
is=find(SREF > 0);
SREF=min(is);
for k=1:length(CORR)
    n(k)=length(CORR(k).CCS);
end
is=find(n > 1);
if min(is) ~= SREF
    SREF=min(is);
end

fprintf('REFERENCE STATION FOR P WAVES IS %s \n',NS(PREF).station);
fprintf('REFERENCE STATION FOR S WAVES IS %s \n',NS(SREF).station);

%%
TTP=(PTEMP-PTEMP(PREF))./SEC; % TEMPI DI ARRIVO P DEL TEMPLATE RELATIVI RISPETTO ALLA REF
TTS=(STEMP-STEMP(SREF))./SEC; % TEMPI DI ARRIVO S DEL TEMPLATE RELATIVI RISPETTO ALLA REF

for k=1:length(NS)
    CORR(k).PCHECK=CORR(k).PT-(TTP(k)*SEC);  %% toglie a tutti i match la travel time dell'evento per le P
    CORR(k).SCHECK=CORR(k).ST-(TTS(k)*SEC);  %% toglie a tutti i match la travel time dell'evento per le S
end

%% ALT CODE FOR EXTRACT EVENTS IMPROVED SPEED
%% Pwaves

PLAG=0.3*SEC; %% search this timelag for P and S
SLAG=0.4*SEC;

fprintf('Collecting P picks\n')
clear POUT
clear ia ib
A=CORR(PREF).PCHECK;
for k=1:length(NS)
    B=CORR(k).PCHECK;
    %[LIA, LIB]=ismembertol(A(1:end),B(1:end),PLAG/(max(abs([A(1:end); B(1:end)]))));
    [LIA, LIB]=ismembertol(A,B,PLAG/(max(abs([A; B]))));
    ia(:,k)=(LIA);
    ib(:,k)=(LIB);
end
%%
S=sum(ia');
ii=find(S>=Args.NP); %% CERCA TUTTI QUELLI CHE HANNO ALMENO NP FASI P
%%
for NEV=1:length(ii)
    n=0;
    TEMP=[];
    for nsta=1:length(NS)
        if (ib(ii(NEV),nsta))~=0
            n=n+1;
            TEMP.sta(n)={NS(nsta).station};    %% real station
            TEMP.Ptim(n)=CORR(nsta).PT(ib(ii(NEV),nsta));  %% real matched P time
            TEMP.Punc(n)=NS(nsta).Punc;
            TEMP.Pcc(n) =CORR(nsta).CCP(ib(ii(NEV),nsta)); %% related cc
            TEMP.Pamp(n)=CORR(nsta).AP(ib(ii(NEV),nsta));
            TEMP.TPamp(n)=max(NS(nsta).ptemplate);
            %% DEBUG
%            TEMP
%            pause
            
        end
    end
    POUT(NEV)=TEMP;
end

if isempty(NEV)==1; fprintf('No P-events found\n');if Args.GARBAGE==1; fprintf('purging working directory...\n');delete('WORK/*');end ;clear all; return; end  %% SE NON TROVA P ESCE
%% S waves
fprintf('Collecting S picks\n')
%SLAG=0.3*SEC;
NEV=0;
clear SOUT
clear ia ib
A=CORR(SREF).SCHECK;

% keyboard
%%
for k=1:length(NS)
    B=CORR(k).SCHECK;
    [LIA, LIB]=ismembertol(A(1:end),B(1:end),SLAG/(max(abs([A(1:end); B(1:end)]))));
    ia(:,k)=(LIA);
    ib(:,k)=(LIB);
end
%%
S=sum(ia');
ii=find(S>=Args.NS); %% CERCA TUTTI QUELLI CHE HANNO ALMENO NP FASI S

for NEV=1:length(ii)
    n=0;
    TEMP=[];
    for nsta=1:length(NS)
        if (ib(ii(NEV),nsta))~=0
            n=n+1;
            TEMP.sta(n)={NS(nsta).station};    %% real station
            TEMP.Stim(n)=CORR(nsta).ST(ib(ii(NEV),nsta));  %% real matched S time
            TEMP.Sunc(n)=NS(nsta).Sunc;
            TEMP.Scc(n) =CORR(nsta).CCS(ib(ii(NEV),nsta)); %% related cc
            TEMP.Samp(n)=CORR(nsta).AS(ib(ii(NEV),nsta));
            TEMP.TSamp(n)=max(NS(nsta).stemplate);
        end
    end
    SOUT(NEV)=TEMP;
end
if isempty(NEV)==1; fprintf('No S-events found\n');if Args.GARBAGE==1; fprintf('purging working directory...\n');delete('WORK/*');end; clear all; return; end

%%
%% MIX P and S matches according to template travel times  ||
fprintf('Combining %8.0d P-events & %8.0d S-events \n',length(POUT),length(SOUT))

%% alt code
SPTIMES=STEMP-PTEMP;  % S-P del template
for k=1:length(NS)
    Ts(k)={NS(k).station};
end
SRCH=0.2/86400;
%keyboard
% INIZIALIZZA SLAVE
SL.Psta=  [];
SL.Ptim=  [];
SL.Punc=  [];
SL.Pcc=   [];
SL.Pamp=  [];
SL.TPamp= [];
SL.Pwvf=  [];
SL.TPwvf= [];
SL.TParr= [];
SL.Ssta=  [];
SL.Stim=  [];
SL.Sunc=  [];
SL.Scc=   [];
SL.Samp=  [];
SL.TSamp= [];
SL.Swvf=  [];
SL.TSwvf= [];
SL.TSarr= [];
SL.TemplateOTime=[];
SL.HEADER=X.HEADER;
SLAVE=repmat(SL,1,5000); %% maximum number of matches per day hardwired!!!
clear SL

for uu=1:length(SOUT)
    Spipo(uu)=SOUT(uu).Stim(1);
end
% keyboard
% Offset is the PREWIN defined at the begin of the program added 17/02/2021
offset=Args.PREWIN/86400;
% keyboard
L=0;
for j=1:length(POUT)
    Pt=POUT(j).Ptim;
    Ps=POUT(j).sta;
    [i, ~]=ismember(Ts,Ps); %% individua fra tutte solo le staz del possibile match
    STEO=Pt+SPTIMES(i);     %% aggiunge alle P del match S-P del template

    %% aggiunto il 4/2/2021 limita la ricerca delle S all'intervallo P -> Steorica+10 secondi
    %MAXS=min(STEO)+(30/86400);
    MAXS=min(STEO)+(2/86400);
    % if MAXS > max(Spipo); MAXS=max(Spipo); end
    %keyboard
    yy=find(Spipo > Pt(PREF) & Spipo < MAXS); % indici
    SOUTTEMP=SOUT(yy);                                   % SOUTTEMP contiene solo le S degli indici
    %%
    % if yy==3; keyboard; end
    for k=1:length(SOUTTEMP)
        St=SOUTTEMP(k).Stim;
        Ss=SOUTTEMP(k).sta;
        NumbS=ismember(Ps,Ss);  %% INDICI DELLE STAZIONI COINCIDENTI
        IS=ismembertol(STEO(NumbS),St,SRCH/(max(abs([STEO(NumbS) St])))); % CERCA GLI ARRIVI S INTORNO ALL'ARRIVO DEL TEMPLATE
        iiS=sum(IS); %% number of STEO ~= St
        %if exist('IS')==1
         if iiS > Args.NS    %% was set to 0 change 28/02/2023
            L=L+1; % keyboard
            SLAVE(L).Psta= POUT(j).sta;
            SLAVE(L).Ptim= POUT(j).Ptim + offset;
            SLAVE(L).Punc= POUT(j).Punc;
            SLAVE(L).Pcc=  POUT(j).Pcc;
            SLAVE(L).Pamp= POUT(j).Pamp;
            SLAVE(L).TPamp= POUT(j).TPamp;            
            SLAVE(L).Ssta= SOUTTEMP(k).sta;
            SLAVE(L).Stim= SOUTTEMP(k).Stim + offset;
            SLAVE(L).Sunc= SOUTTEMP(k).Sunc;
            SLAVE(L).Scc=  SOUTTEMP(k).Scc;
            SLAVE(L).Samp= SOUTTEMP(k).Samp;
            SLAVE(L).TSamp=SOUTTEMP(k).TSamp;
            SLAVE(L).TemplateOTime=NS(1).MOTime;
%             %% DEBUG
%             SLAVE(L)
%             pause
        end
    end
end
SLAVE=SLAVE(1:L); %% RESIZE SLAVE

fprintf('===============================\n');
fprintf('# total matches found= %4.0f   |\n',length(SLAVE));
fprintf('===============================\n');

%% FILL SLAVES WITH WAVEFORMS
%% STARTING PRE SEC FROM P TIME 
%% ADDED 11/02/2023
PRE=Args.PREWIN+1; PPOS=Args.PLWIN+1; SPOS=Args.SLWIN+1;
for k=1:length(FULLDAYIDx)    
    STAFULLDAY{k}=deblank(FULLDAYIDx{k}.station);
end
for k=1:length(NS)    
    STANS{k}=deblank(NS(k).station);
end

for k=1:L
    for ind=1:numel(SLAVE(k).Ptim);
        STA =SLAVE(k).Psta(ind);
        ix=find(contains(STAFULLDAY,STA)); %% FIND THE  RIGHT STATION IN FULLDAY
        [minDistance, indexOfMin] = min(abs(FULLDAYIDx{ix}.matlabTimeVector-SLAVE(k).Ptim(ind)));
        SMP=FULLDAYIDx{ix}.sampleRate;
        SLAVE(k).Pwvf(ind).time=FULLDAYIDx{ix}.matlabTimeVector(indexOfMin-(round(PRE*SMP)):indexOfMin+(round(PPOS*SMP)));
        SLAVE(k).Pwvf(ind).data=FULLDAYIDx{ix}.datafilt(indexOfMin-(round(PRE*SMP)):indexOfMin+(round(PPOS*SMP)));
        
        ins=find(contains(STANS,STA)); %% FIND THE  RIGHT STATION IN NS
        [minDistance, indexOfMin] = min(abs(NS(ins).matlabTimeVector-NS(ins).PTIME));
        SLAVE(k).TPwvf(ind).time=NS(ins).matlabTimeVector(indexOfMin-(round(PRE*SMP)):indexOfMin+(round(PPOS*SMP)));
        SLAVE(k).TPwvf(ind).data=NS(ins).datafilt(indexOfMin-(round(PRE*SMP)):indexOfMin+(round(PPOS*SMP)));
        SLAVE(k).TParr(ind)=NS(ins).PTIME;
    end
    for ind=1:numel(SLAVE(k).Stim);
        STA =SLAVE(k).Ssta(ind);
        ix=find(contains(STAFULLDAY,STA)); %% FIND THE  RIGHT STATION IN FULLDAY
        [minDistance, indexOfMin] = min(abs(FULLDAYIDx{ix}.matlabTimeVector-SLAVE(k).Stim(ind)));
        SMP=FULLDAYIDx{ix}.sampleRate;
        SLAVE(k).Swvf(ind).time=FULLDAYIDx{ix}.matlabTimeVector(indexOfMin-(round(PRE*SMP)):indexOfMin+(round(SPOS*SMP)));
        SLAVE(k).Swvf(ind).data=FULLDAYIDx{ix}.datafilt(indexOfMin-(round(PRE*SMP)):indexOfMin+(round(SPOS*SMP)));
        
        ins=find(contains(STANS,STA)); %% FIND THE  RIGHT STATION IN NS
        [minDistance, indexOfMin] = min(abs(NS(ins).matlabTimeVector-NS(ins).STIME));
        SLAVE(k).TSwvf(ind).time=NS(ins).matlabTimeVector(indexOfMin-(round(PRE*SMP)):indexOfMin+(round(SPOS*SMP)));
        SLAVE(k).TSwvf(ind).data=NS(ins).datafilt(indexOfMin-(round(PRE*SMP)):indexOfMin+(round(SPOS*SMP)));
        SLAVE(k).TSarr(ind)=NS(ins).STIME;
    end        
end

ET=etime(clock,T0);
fprintf('Total time used for collect %4.0f Events: %8.2f sec \n',length(SLAVE),ET)

save(['OUT/' OTime '.' Args.DAYDIR '.mat'],'SLAVE');

% for i=1:length(SLAVE); fprintf('%d  %d \n',length(SLAVE(i).Pcc)-length(SLAVE(i).Punc),length(SLAVE(i).Scc)-length(SLAVE(i).Sunc)); end

%% CLEAR OUT WORKFILES
if Args.GARBAGE==1
    fprintf('purging working directory...\n')
    !rm WORK/*
end

%parfevalOnAll(@clearvars, 0)
clear SLAVE
% clear all
%% END OF MAIN

