%% Prepara il file dt.cc leggendo in ingresso event.sel e dt.cc di PH2DT
%% Le f.o. sono sorganizzate in dir ordinate per ID
%% 


clear all
close all

% write or not data
ws=0;

% READ LOCAL DATA or NOT
rl=0;
if rl==1;
    fprintf('Reading data from disk\n')
end

% cc threshold
THR=0.7;

% sac file length
PRE= 3/86400;
DUR=25/86400;

FMT='yyyy-mm-ddTHH:MM:SS.FFF';
FOUT='YY.XXXX.ZZZ'; 
filo=fopen('dt.cc','a');  %% SET FILEOUT NAME


%% MAIN
% Read in catalog and IDS
CT=textread('./dt.ct','%s','delimiter','\n');
idx=find(contains(CT,'#'));


% THESE EVENTS ARE SELECTED BY PH2DT
% Load event.sel
EV=readtable('../PH2DT/event.sel','FileType','text');
YMD=num2str(EV.Var1);
HMS=num2str(EV.Var2);

for k=1:length(YMD)
    YMD(k,:)=strrep(YMD(k,:),' ','0');
    HMS(k,:)=strrep(HMS(k,:),' ','0');
    OTIME(k)=datenum([YMD(k,:) HMS(k,1:6) '.' HMS(k,7:8) '0'],'yyyymmddHHMMSS.FFF');
end




for k=1:numel(idx)-1    
    BLK=CT(idx(k):idx(k+1)-1);
    [junk,ID1,ID2]=strread(char(BLK(1)),'%s %d %d');
    fprintf(filo,'#    %4d   %4d     0\n',ID1,ID2);
    
    fprintf('#    %4d   %4d    %d/%d\n',ID1,ID2,k,numel(idx));

    if ws==1
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('./',sprintf('%d',ID1));
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('./',sprintf('%d',ID2));
    end
    
    i1=find(EV.Var10==ID1);
    i2=find(EV.Var10==ID2);
    OT1=OTIME(i1);
    OT2=OTIME(i2);
    
    TIMEO1pre=datestr(OT1-PRE,FMT);
    TIMEO1pos=datestr(OT1+DUR,FMT);
    TIMEO2pre=datestr(OT2-PRE,FMT);
    TIMEO2pos=datestr(OT2+DUR,FMT);
    
    % keyboard
    BLK=BLK(2:end);
    S=[];
    for n=1:length(BLK)    
        [STA,P1,P2,W,PHS]=strread(char(BLK(n)),'%s %f %f %f %s');
        if isempty (str2num(char(STA)))==0; STA={['T' char(STA)]}; end
%        i=find(contains(BLK,STA));
        
        if isfield(S,char(STA))==0
            S.(char(STA))=[];
            S.(char(STA)).ID1=ID1;
            S.(char(STA)).ID2=ID2;
            S.(char(STA)).OT1=OT1;
            S.(char(STA)).OT2=OT2;
            S.(char(STA)).P1=[];
            S.(char(STA)).P2=[];
            S.(char(STA)).S1=[];
            S.(char(STA)).S2=[];
            S.(char(STA)).Pw=[];
            S.(char(STA)).Sw=[];
        end
        if strcmp(PHS,'P')
            S.(char(STA)).P1=P1;
            S.(char(STA)).P2=P2;
            S.(char(STA)).Pw=W;
        else
            S.(char(STA)).S1=P1;
            S.(char(STA)).S2=P2;
            S.(char(STA)).Sw=W;
        end
    end
            
    NSTA=fieldnames(S);
    
   
    % S is a structure containing each information        
    for j=1:numel(NSTA)
        STA=char(NSTA(j));
        
        if rl==0 % DOWNLOAD DATA
            %% ID1
            P1=S.(STA).P1;
            S1=S.(STA).S1;
            [MTIME,DATA,OUT1,WVF]=get_waveforms('INGV',TIMEO1pre ,TIMEO1pos ,'*',STA,'*','HHZ');
            if isempty(OUT1)==1
                [MTIME,DATA,OUT1,WVF]=get_waveforms('INGV',TIMEO1pre ,TIMEO1pos ,'*',STA,'*','EHZ');
            end
            if isempty(OUT1)==1
                [MTIME,DATA,OUT1,WVF]=get_waveforms('INGV',TIMEO1pre ,TIMEO1pos ,'*',STA,'*','CHZ');
            end
            
            if isempty(OUT1)==1 
               NODATA1=1;
            else
                NODATA1=0;
            end
            
            if NODATA1==0
            NET=char(WVF.traces{1}.stats.network);
            STA=char(WVF.traces{1}.stats.station);
            CMP=char(WVF.traces{1}.stats.channel);
            FOUT=[NET '.' STA '.' CMP];
            if P1 > 0
                OUT1=ch(OUT1,'A',P1+(PRE*86400));
            end
            if S1 > 0
                OUT1=ch(OUT1,'T0',S1+(PRE*86400));
            end
            OUT1=ch(OUT1,'O',(PRE*86400));
            OUT1=ch(OUT1,'KNETWK',NET);
            if ws==1
                wsac([sprintf('%d/',ID1) FOUT],OUT1);
            end
            end
            
            %% ID2
            P2=S.(STA).P2;
            S2=S.(STA).S2;
            [MTIME,DATA,OUT2,WVF]=get_waveforms('INGV',TIMEO2pre ,TIMEO2pos ,'*',STA,'*','HHZ');
            if isempty(OUT2)==1
                [MTIME,DATA,OUT2,WVF]=get_waveforms('INGV',TIMEO2pre ,TIMEO2pos ,'*',STA,'*','EHZ');
            end
            if isempty(OUT2)==1
                [MTIME,DATA,OUT2,WVF]=get_waveforms('INGV',TIMEO2pre ,TIMEO2pos ,'*',STA,'*','CHZ');
            end
            
            if isempty(OUT2)==1 
               NODATA2=1;
            else
                NODATA2=0;
            end
            
            if NODATA2==0
            NET=char(WVF.traces{1}.stats.network);
            STA=char(WVF.traces{1}.stats.station);
            CMP=char(WVF.traces{1}.stats.channel);
            FOUT=[NET '.' STA '.' CMP];
            if P2 > 0
                OUT2=ch(OUT2,'A',P2+(PRE*86400));
            end
            if S2 > 0
                OUT2=ch(OUT2,'T0',S2+(PRE*86400));
            end
            OUT2=ch(OUT2,'O',(PRE*86400));
            OUT2=ch(OUT2,'KNETWK',NET);
            if ws==1
                wsac([sprintf('%d/',ID2) FOUT],OUT2);
            end
            end
            
        else
            if NODATA1==0 && NODATA2==0
            D1=dir([sprintf('%d/*',ID1) STA '.*']);
            OUT1=rsac([D1.folder '/' D1.name]);
            
            D2=dir([sprintf('%d/*',ID2) STA '.*']);
            OUT2=rsac([D2.folder '/' D2.name]);
            
            P1=S.(STA).P1;
            S1=S.(STA).S1;
            P2=S.(STA).P2;
            S2=S.(STA).S2;
            end                  
        end
        
        if NODATA1==0 & NODATA2==0
        if P1 > 0
            [FINECC,FINESHIFT,PHS]=do_cc(OUT1,OUT2,S.(STA).P1,S.(STA).P2,'P',[2 15]);
            if FINECC > THR;
                %STA=STA(2:end);
                fprintf(filo,'%5s % 6.4f  %4.2f  %s\n',STA,FINESHIFT,FINECC,'P');
            end
        end
        if S1 > 0                    
            [FINECC,FINESHIFT,PHS]=do_cc(OUT1,OUT2,S.(STA).S1,S.(STA).S2,'S',[2 15]);
            if FINECC > THR;
                %STA=STA(2:end);
                fprintf(filo,'%5s % 6.4f  %4.2f  %s\n',STA,FINESHIFT,FINECC,'S');
            end
        end   
        end
        
    end        
end              
fclose(filo);

%% FUNCTIONS

function [FINECC,FINESHIFT,PHS]=do_cc(sac1,sac2,t1,t2,PHS,filt)


FIG=0;


if strcmp(PHS,'P')
    PRE=0.2;
    POS=1.0;
    MAXLAG=.8;
else
    PRE=0.5;
    POS=1.5;
    MAXLAG=1.0;
end


ML=round(MAXLAG*100);

sac1=condiziona(sac1,filt);
sac2=condiziona(sac2,filt);

if FIG==1
figure;
p1(sac1);
figure;
p1(sac2);
end

O1=lh(sac1,'O');
O2=lh(sac2,'O');

%% RESIZE P LENGTH
if strcmp(PHS,'P')
SP=lh(sac1,'T0')-lh(sac1,'A');
if SP > 0 && SP < POS
    % resize winlen        
    POS=2*SP/3;
    MAXLAG=PRE+POS;
    fprintf('Resizing win len to %4.2f s\n',POS)
end
    
end

%%


i1=find(sac1(:,1) > t1-PRE-lh(sac1,'B')+O1 & sac1(:,1) < t1+POS-lh(sac1,'B')+O1);
i2=find(sac2(:,1) > t2-PRE-lh(sac2,'B')+O2 & sac2(:,1) < t2+POS-lh(sac2,'B')+O2);

L1=numel(i1);
L2=numel(i2);

if L1 ~= L2
    if L1 > L2;
        i1=i1(1:L2);
    else
        i2=i2(1:L1);
    end
end

   
if FIG==1
    figure;
    subplot(211)
    plot(normize(sac1(i1,2)));
    hold on
    plot(normize(sac2(i2,2)));
end

[C, LAG]=xcorr(sac1(i1,2),sac2(i2,2),ML,'normalized'); 
if FIG==1; subplot(212); plot(LAG,C); end
[FINESHIFT,FINECC]=subsample(LAG,C);
FINESHIFT=-FINESHIFT*lh(sac1,'DELTA'); %% <<-- corretto il 12 febbraio 2023 su SIENA
% fprintf('%6.4f  %4.2f  %s\n',FINESHIFT,FINECC,PHS);

end

function [FINESHIFT,FINECC]=subsample(LAG,C)


FIT=1;  % 0=cubic, 1=spline

RSMP=0.001;

i=find(C==max(C)); 
if i > 2 & i < numel(C)-2
    FINECC=C(i);
    ILAG=LAG(i);
    PORT=C(i-2:i+2);
    if FIT==0
        P = polyfit(ILAG-2:ILAG+2,PORT',2);
        X=ILAG-2:RSMP:ILAG+2;
        Y=((X.^2)*P(1))+(X*P(2))+P(3);
    else
        X =ILAG-2:ILAG+2;
        XQ=ILAG-2:RSMP:ILAG+2;
        Y=spline(X,PORT',XQ);
    end
    NSUB=((numel(X)-1)/2)-find(Y==max(Y));
    FINESHIFT=ILAG - (NSUB*RSMP);
else
    FINESHIFT=0;
    FINECC   =0;
end
    
    
end













