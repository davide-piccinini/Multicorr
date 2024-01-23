function []=selectcorr(configfile)

% USAGE []=selectcorr(configfile);
%
% Makes some pretty figures and
% SELECT CORRELATION OUTPUT according to
% CC threshold and NP / NS
% then write out NLL obs files ready to be relocated using
% NLLoc (Lomax)
% To be used after mcorr
%
% 17/10/2022 - Davide Piccinini
% + MAD threshold type was implemented
%
% threshold is calculated as CFG.MAD * (mad(Pcc))
%
% 17/10/2022
% Amplitude based Magnitude estimation implemented
% M2=M1+log(A2)-log(A1)


%% Read Configuration File
if nargin < 1
    % DEFAULT PARAMETERS
    CFG.MAD     = 10;
    CFG.CCTHR   = 0.73;
    CFG.SWITCH  = 0;
    CFG.NP      = 5;
    CFG.NS      = 4;
    CFG.MATDIR  = 'OUT/';
    CFG.NLLOUT  = 'NLL_OBS/';
    CFG.DT      = 5;
    CFG.FIG     = 1;
else
    CFG=readcfg(configfile);  %% Now use the readcfg.m function (added 5 oct 2022)
end

NF=fieldnames(CFG);
fprintf('Running using following parameters:\n')
for k=1:length(NF);
    VAR=getfield(CFG,char(NF(k)));
    fprintf('%s = %s\n',char(NF(k)),string(VAR))
end


%% Check if NLL DIRECTORY exist
EXI=exist(CFG.NLLOUT);
if EXI==0
    mkdir(CFG.NLLOUT);
end

%keyboard
%% Work
D=dir([CFG.MATDIR '*.mat']);


L=0;
THR=CFG.CCTHR;

for k=1:length(D)
    S=load([CFG.MATDIR D(k).name]);
    N=length(S.SLAVE);
    if N > 0; % SLAVE is not empty
        for j=1:N
            NUMP=length(S.SLAVE(j).Pcc);
            NUMS=length(S.SLAVE(j).Scc);
            MP= mean(S.SLAVE(j).Pcc);
            MS= mean(S.SLAVE(j).Scc);
            MTP=(S.SLAVE(j).Ptim(1));
            MTS=(S.SLAVE(j).Stim(1));
            
            MADP=CFG.MAD*mad(S.SLAVE(j).Pcc,1);
            MADS=CFG.MAD*mad(S.SLAVE(j).Scc,1);
            
            if CFG.SWITCH==0
                LOGIC=MP >= THR && MS >= THR && NUMP >= CFG.NP && NUMS >= CFG.NS;
            else
                LOGIC=CFG.CCTHR >=MADP && CFG.CCTHR>=MADS && NUMP >= CFG.NP && NUMS >= CFG.NS;
            end
            
            if LOGIC==1
                try
                    L=L+1;
                    SLA(L).TemplateOTime=S.SLAVE(j).TemplateOTime;
                    SLA(L).Psta= S.SLAVE(j).Psta;
                    SLA(L).Ptim= S.SLAVE(j).Ptim;
                    SLA(L).Pcc=  S.SLAVE(j).Pcc;
                    SLA(L).Ssta= S.SLAVE(j).Ssta;
                    SLA(L).Stim= S.SLAVE(j).Stim;
                    SLA(L).Scc=  S.SLAVE(j).Scc;
                    SLA(L).MTP=MTP;
                    SLA(L).MTS=MTS;
                    SLA(L).MP=MP;
                    SLA(L).MS=MS;
                    SLA(L).Pamp=  S.SLAVE(j).Pamp;
                    SLA(L).TPamp= S.SLAVE(j).TPamp;
                    SLA(L).Samp=  S.SLAVE(j).Samp;
                    SLA(L).TSamp= S.SLAVE(j).TSamp;
                    SLA(L).HEADER=S.SLAVE(j).HEADER;
                    SLA(L).MAGP   =mean(SLA(L).HEADER{6}+log10(SLA(L).Pamp)-log10(SLA(L).TPamp));
                    SLA(L).MAGS   =mean(SLA(L).HEADER{6}+log10(SLA(L).Samp)-log10(SLA(L).TSamp));
                    SLA(L).MAG=mean([SLA(L).MAGP SLA(L).MAGS]);
                    SLA(L).MAGERR=std([SLA(L).MAGP SLA(L).MAGS]);
                catch
                    fprintf("Houston we had a problem here...\n")
                    keyboard
                    %return
                end
                
            end
        end
    end
end

fprintf('Total # of events with CC > %3.2f & NP >= %d & NS >= %d : %d\n',CFG.CCTHR,CFG.NP,CFG.NS,length(SLA))

T=struct2table(SLA);
%keyboard
% S=sortrows(T,8);
%%%%%%%%%%%%%
TT=unique(T(:,'TemplateOTime'));
for k=1:height(T)
    PTIME(k)=T{k,'Ptim'}{1}(1);
end

for k=1:height(TT)
    OT=TT{k,'TemplateOTime'};
    i=find(PTIME > OT & PTIME < OT+(10/86400));
    if isempty(i)==0
        datestr(OT);
        datestr(PTIME(i));
        IDK=find((T{i,'TemplateOTime'})==OT);
        if length(IDK) > 0
            MATCH(k)=i(IDK(1));
            %pause
        end
    end
end
% ALCUNI TEMPLATES POTREBBERO NON AVERE MATCHES (STRANO)
MATCH(find(MATCH==0))=[]; %% TOGLIE QUELLI CHE NON HANNO NEANCHE UN MATCH

%keyboard
TEMPLATES=table2struct(T(MATCH,:)); % QUESTI SONO I TEMPLATES
T(MATCH,:)=[];                      % E LI TOLGO DAL TOTALE
PTIME(MATCH)=[];
%%%%%%%%%%%%%
%    TIMEP=T{:,8};
[PTIME idx]=sort(PTIME); %% rearrange Ptime and table according to ascending Ptime
T=T(idx,:);

TOL=CFG.DT/86400;
[C,ID,IJ] = uniquetol(PTIME,TOL/max(abs(PTIME)));

NMAX=max(IJ);
TOERASE=[];

for j=1:NMAX
    NSTAP=[];NSTAS=[];CCP=[];CCS=[];
    II=find(IJ==j);
    for k=1:length(II)
        NSTAP(k)=length(T{II(k),'Psta'}{:});
        NSTAS(k)=length(T{II(k),'Ssta'}{:});
        CCP  (k)=(T{II(k),'MP'});
        CCS  (k)=(T{II(k),'MS'});
    end
    datestr(C(j));
    Z=[NSTAP ; NSTAS; CCP; CCS];
    SCORE=sum(Z);
    [w v]=max(SCORE);
    CHOICE=II(v);
    
    OUT=setdiff(II,CHOICE);
    LO=length(OUT);
    LE=length(TOERASE);
    TOERASE(LE+1:LE+length(OUT))=OUT ;
    % pause
end

[r c]=size(T);
if exist('TOERASE','var') == 1
    T([TOERASE],:)=[]; % TOGLIE GLI INDICI CHE HA SCARTATO
    fprintf('%4.0f double events removed\n',length(TOERASE))
else
    fprintf('No doubles found, %4d events keeped\n',r);
end
SLAVE=vertcat(TEMPLATES,table2struct(T));

[r,c]=size(SLAVE);
% fprintf('Total slaves collected: %d\n',r)
clear SLA TOERASE II

fprintf('TOTAL NUMBER OF SLAVES COLLECTED: %4d\n',r)
fprintf('NUMBER OF TEMPLATES FOUND       : %4d\n', length(TEMPLATES))
fprintf('NUMBER OF NEW EVENTS FOUND      : %4d\n', height(T))

save('RESULTS.mat','SLAVE')
% keyboard



%%
% PLOT SERIES
for k=1:length(SLAVE); OTIME(k)=SLAVE(k).TemplateOTime;TIMEPP(k)=SLAVE(k).Ptim(1);TIMES(k)=SLAVE(k).Stim(1);CCP(k)=SLAVE(k).MP; end

U=unique(OTIME);  %% U IS THE TEMPLATE OTIME
clear TU TS
clear CC

% for k=1:length(U);
%     i=find(OTIME==U(k));
%     TU(k,2:length(i)+1)=NaN;
%     TU(k,1)=TIMEPP(k);
%     TS(k,2:length(i)+1)=NaN;
%     TS(k,1)=TIMES(k);
%     CC(k,2:length(i)+1)=NaN;
%     CC(k,1)=CCP(k);
%
%     for j=1:length(i)
%         TU(k,j+1)=TIMEPP(i(j));
%         TS(k,j+1)=TIMES(i(j));
%         CC(k,j+1)=CCP(i(j));
%     end
%     i=[];
% end

%%
TU=[];
TS=[];
CC=[];
for k=1:length(U)
    ix=find(OTIME==U(k));
    Lu=length(TU);
    Li=length(ix);
    TU(k,Lu+1:Lu+Li)=TIMEPP(ix);
    CC(k,Lu+1:Lu+Li)=CCP(ix);
    
end


[r c]=size(TU); %% THIS IS THE SCATTERPLOT VARIABLE

% EVtime(1)=datenum(2022,5,3,15,50,49);EVM(1)=3.7;
% EVtime(2)=datenum(2022,5,3,20,14,20);EVM(2)=3.5;
% EVtime(3)=datenum(2022,5,10,3,51,17);EVM(3)=3.5;
% EVtime(4)=datenum(2022,5,12,21,12,3);EVM(4)=3.7;

% keyboard
%%
if CFG.FIG==1
    figure;
    subplot(211)
    for k=1:r
        ii=find(TU(k,:) > 0);
        X=k.*ones(size(TU(k,ii)));
        % scatter(TU(k,ii(2:end)),X(2:end),30,CC(k,ii(2:end)),'s','filled'); hold on
        scatter(TU(k,ii),X,30,CC(k,ii),'s','filled'); hold on
    end
    datetick('x','mmdd','keeplimits')
    box on
    grid on
    ylabel('Template ID')
    colorbar('East')
    
    % hold on
    % scatter(EVtime,ones(1,4)*r,EVM.^4,'pk','MarkerFaceColor','y')
    
    
    %%
    W=zeros(1,c);
    W(1:r)=1;
    W=[W ; zeros(size(W)); zeros(size(W))];W=W';
    %%
    [TIMEPP idx]=sort(TIMEPP);
    W=W(idx,:);
    %%
    subplot(212);
    scatter(TIMEPP,1:length(TIMEPP),20,[W],'filled')
    %plot(TIMEPP,1:length(TIMEPP),'sk','MarkerFaceColor','r','MarkerSize',4);
    datetick('x','mm/dd','keeplimits');grid on
    xlabel('Time ')
    ylabel('Number of events')
    subplot(212); title(sprintf('%s - %s - NP=%d NS=%d CC=%4.2f',datestr(floor(TIMEPP(1)),'yyyymmddHHMM'),datestr(ceil(TIMEPP(end)),'yyyymmddHHMM'),CFG.NP,CFG.NS,CFG.CCTHR))
    
    % hold on
    % scatter(EVtime,ones(1,4)*100,EVM.^4,'pk','MarkerFaceColor','y')
end

%% % WRITE OUT NNLINLOC PHS FILE
fprintf('Removing old nlloc observation files...\n')
delete [CFG.NLLOUT '*']
fprintf('Writing new nlloc observation files...\n')
STHR=sprintf('%3.2f',THR);

for k=1:length(SLAVE)
    SLAVE(k).Psta;
    NAME=datestr(SLAVE(k).Ptim(1),'yyyymmdd_HHMMSS.FFF');
    fout=fopen([CFG.NLLOUT NAME '.' STHR '.xml'],'w');
    for j=1:length(SLAVE(k).Psta);
        PSTA=SLAVE(k).Psta(j);
        PCH ='Z';
        PHS ='P';
        TIME=datestr(SLAVE(k).Ptim(j),'yyyymmdd HHMM SS.FFF');
        WW  =(1-SLAVE(k).Pcc(j))/3;
        fprintf(fout,'%-6s ?    %s    ? %s      ? %s0 GAU %9.2e -1.00e+00 -1.00e+00 -1.00e+00\n',...
            char(PSTA),PCH,PHS,TIME,WW);
    end
    for j=1:length(SLAVE(k).Ssta);
        SSTA=SLAVE(k).Ssta(j);
        SCH ='Z';
        PHS ='S';
        TIME=datestr(SLAVE(k).Stim(j),'yyyymmdd HHMM SS.FFF');
        WW  =(1-SLAVE(k).Scc(j))/3;
        fprintf(fout,'%-6s ?    %s    ? %s      ? %s0 GAU %9.2e -1.00e+00 -1.00e+00 -1.00e+00\n',...
            char(SSTA),SCH,PHS,TIME,WW);
    end
    fclose(fout);
end

% %%
% ZIPFILE=sprintf('zip IMPRUNETA_%1dp_%1ds_%4.2f.zip NLL_OBS/*.xml',NP,NS,THR);
% eval(ZIPFILE)


%% print out the N most productive templates
%% this could help to download waveforms for a couple template/station
FMT='yyyy-mm-ddTHH:MM:SS.FFF';

N=3;
for k=1:r; % size(TU);
    index(k)=length(find(TU(k,:) > 0 ));
end

[inds i]=sort(index,'descend');
TUsort=TU(i,:);
Usort =U(i);

fprintf('These are the most %2.0f productive templates:\n',N)
fprintf(' \n',N)

for k=1:N
    idx=find(TUsort(k,:) > 0);
    mint=min(TUsort(k,idx));
    maxt=max(TUsort(k,idx));
    
    fprintf('Template %2.0f %3.0f %s    %s -> %s \n',k,length(idx),datestr(Usort(k),FMT),datestr(mint,FMT),datestr(maxt,FMT));
end

% keyboard
%% do you want to download all the matched waveforms for selected sta/template??
%% this is an example for the most productive
% ID=4;
% FMT='yyyy-mm-ddTHH:MM:SS.FFF';
%
% STA='HLNI';
% NET='IV';
% CHA='HHZ';
% PRE=10;
% POS=30;
%
% NM=sprintf('%s',datestr(TU(ID),FMT))
% [SUCCESS,MESSAGE,MESSAGEID] = mkdir(NM)
%
% idx=find(TU(ID,:) > 0);
% for k=1:length(idx);
%     TP=TU(ID,idx(k));
%     TINI=datestr(TP-(PRE/86400),FMT);
%     TFIN=datestr(TP+(POS/86400),FMT);
%     [TIME D O W]=get_waveforms('INGV',TINI,TFIN,'*',STA,'',CHA);
%     wsac([NM '/' datestr(TP,FMT) '.' STA '.' CHA],O);
% end


%
% D=dir([NM '/*HHZ']);
%
% for k=1:15
%     A=rsac([D(k).folder '/' D(k).name]);
%     if k==1;
%         STACK=A(:,2);
%     else
%         STACK=STACK+A(:,2);
%     end
% end








    
    %% END OF MAIN
    

function CFG=readconfig(configfile);

f=fopen(configfile); A=textscan(f,'%s','Delimiter','\n'); fclose(f);
A=char(A{1});
CFG.CCTHR=str2num(A(2,:));
CFG.NP=str2num(A(4,:));
CFG.NS=str2num(A(6,:));
CFG.MATDIR=deblank(A(8,:));
CFG.NLLOUT=deblank(A(10,:));
CFG.DT=str2num(A(12,:));
CFG.FIG=str2num(A(14,:));
