% function correlate
clear all
% cc threshold
THR=0.6;

FMT='yyyy-mm-ddTHH:MM:SS.FFF';
FOUT='YY.XXXX.ZZZ';
filo=fopen('dt.cc','a');  %% SET FILEOUT NAME

%%
load /homes/piccinini/MULTICORR/SIENA/RESULTS.mat

for k=1:numel(SLAVE)
    SLTIME(k)=SLAVE(k).Ptim(1);
end


%% MAIN
% Read in catalog and IDS
CT=textread('/homes/piccinini/MULTICORR/SIENA/PH2DT/dt.ct','%s','delimiter','\n');
idx=find(contains(CT,'#'));


% THESE EVENTS ARE SELECTED BY PH2DT
% Load event.sel
EV=readtable('/homes/piccinini/MULTICORR/SIENA/PH2DT/event.sel','FileType','text');
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
    
    
    i1=find(EV.Var10==ID1);
    i2=find(EV.Var10==ID2);
    OT1=OTIME(i1);
    OT2=OTIME(i2);
    
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
    
    ID_SLAVE_1=ID1;
    ID_SLAVE_2=ID2;
    
    
    
    %     TTP1=S.(char(NSTA(1))).P1;
    %     if isempty(TTP1);
    %         TTP1=S.(char(NSTA(1))).S1;
    %     end
    %     TTP1=TTP1/86400;
    %     PAIROT1=OT1+TTP1;
    %     ID_SLAVE_1=find(abs(SLTIME-PAIROT1) < 9/86400);
    %
    %     TTP2=S.(char(NSTA(1))).P2;
    %     if isempty(TTP2);
    %         TTP2=S.(char(NSTA(1))).S2;
    %     end
    %     TTP2=TTP2/86400;
    %     PAIROT2=OT2+TTP2;
    %     ID_SLAVE_2=find(abs(SLTIME-PAIROT2) < 9/86400);
    % if ID1==623 && ID2==624; keyboard; end
    
    for jj=1:numel(NSTA)
        S.(char(NSTA(jj)));
        STA=char(NSTA(jj));
        P1=S.(char(NSTA(jj))).P1;
        if isempty(P1)==0
            ik1=find(contains(SLAVE(ID_SLAVE_1).Psta,STA));
            ik2=find(contains(SLAVE(ID_SLAVE_2).Psta,STA));
            OUT1.data=SLAVE(ID_SLAVE_1).Pwvf(ik1).data;
            OUT2.data=SLAVE(ID_SLAVE_2).Pwvf(ik2).data;
            OUT1.time=SLAVE(ID_SLAVE_1).Pwvf(ik1).time;
            OUT2.time=SLAVE(ID_SLAVE_2).Pwvf(ik2).time;
            OUT1.OT=OT1;
            OUT2.OT=OT2;
            OUT1.TTP=S.(STA).P1;
            OUT2.TTP=S.(STA).P2;
            OUT1.TTS=S.(STA).S1;
            OUT2.TTS=S.(STA).S2;
            [FINECC,FINESHIFT,PHS]=do_cc(OUT1,OUT2,'P');
            if FINECC >= THR;
                fprintf(filo,'%5s % 6.4f  %4.2f  %s\n',STA,FINESHIFT,FINECC,'P');
            end
        end
        
        S1=S.(char(NSTA(jj))).S1;
        if isempty(S1)==0
            ik1=find(contains(SLAVE(ID_SLAVE_1).Ssta,STA));
            ik2=find(contains(SLAVE(ID_SLAVE_2).Ssta,STA));
            OUT1.data=SLAVE(ID_SLAVE_1).Swvf(ik1).data;
            OUT2.data=SLAVE(ID_SLAVE_2).Swvf(ik2).data;
            OUT1.time=SLAVE(ID_SLAVE_1).Swvf(ik1).time;
            OUT2.time=SLAVE(ID_SLAVE_2).Swvf(ik2).time;
            OUT1.OT=OT1;
            OUT2.OT=OT2;
            OUT1.TTP=S.(STA).P1;
            OUT2.TTP=S.(STA).P2;
            OUT1.TTS=S.(STA).S1;
            OUT2.TTS=S.(STA).S2;
            [FINECC,FINESHIFT,PHS]=do_cc(OUT1,OUT2,'S');
            % fprintf('%5s % 6.4f  %4.2f  %s\n',STA,FINESHIFT,FINECC,'S');
            if FINECC >= THR;
                fprintf(filo,'%5s % 6.4f  %4.2f  %s\n',STA,FINESHIFT,FINECC,'S');
            end
        end
    end
    
end

    
    
    
    
    
    
    %% FUNCTIONS

function [FINECC,FINESHIFT,PHS]=do_cc(OUT1,OUT2,PHS)


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

if FIG==1
    if OUT1.TTP > 0 | OUT1.TTS > 0  
    figure
    subplot(211); plot(OUT1.time,OUT1.data); hold on; plot(OUT1.OT,0,'*r'); plot(OUT1.OT+(OUT1.TTP/86400),0,'+k'); datetick('x','MM:SS.FFF')
    subplot(212); plot(OUT2.time,OUT2.data); hold on; plot(OUT2.OT,0,'*r'); plot(OUT2.OT+(OUT2.TTP/86400),0,'+k'); datetick('x','MM:SS.FFF')
    end
end

%% RESIZE P LENGTH
if strcmp(PHS,'P') & OUT1.TTS > 0
    SP=OUT1.TTS-OUT1.TTP;
    if SP > 0 & SP < POS
        % resize winlen
        POS=2*SP/3;
        MAXLAG=PRE+POS;
        fprintf('Resizing win len to %4.2f s\n',POS)
    end
end

%%
if strcmp(PHS,'P')
    t1=OUT1.TTP;
    t2=OUT2.TTP;
else
    t1=OUT1.TTS;
    t2=OUT2.TTS;
end


i1=find(OUT1.time > OUT1.OT+(t1/86400)-(PRE/86400) & OUT1.time < OUT1.OT+(t1/86400)+POS);
i2=find(OUT2.time > OUT2.OT+(t2/86400)-(PRE/86400) & OUT2.time < OUT2.OT+(t2/86400)+POS);

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
    if t1 > 0
    figure;
    subplot(211)
    plot(normize(OUT1.data(i1)));
    hold on
    plot(normize(OUT2.data(i2)));axis tight
    end
end

[C, LAG]=xcorr(OUT1.data(i1),OUT2.data(i2),ML,'normalized');
if FIG==1; subplot(212); plot(LAG,C); end

DELTA=mean(diff(OUT1.time))*86400; SMP=round(1/DELTA);
[FINESHIFT,FINECC]=subsample(LAG,C);
FINESHIFT=-FINESHIFT*DELTA; %% <<-- corretto il 12 febbraio 2023 su SIENA
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


