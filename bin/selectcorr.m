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
%
% Last valid modification Apr 28 2023

% CFG.MULMAG added Jan 3 2024

%% Read Configuration File
if nargin < 1
    % DEFAULT PARAMETERS
    CFG.MAD     = 14;
    CFG.CCTHR   = 0.45;
    CFG.SWITCH  = 0;
    CFG.NP      = 4;
    CFG.NS      = 4;
    CFG.MATDIR  = 'OUT/';
    CFG.NLLOUT  = 'NLL_OBS/';
    CFG.DT      = 10;
    CFG.FIG     = 1;
    CFG.JSON    = 0;
    CFG.MAGFILT = 1;
    CFG.MULMAG  = 1;
else
    CFG=readcfg(configfile);  %% Now use the readcfg.m function (added 5 oct 2022)
end

NF=fieldnames(CFG);
fprintf('Running using following parameters:\n')
for k=1:length(NF);
    VAR=getfield(CFG,char(NF(k)));
    fprintf('%s = %s\n',char(NF(k)),string(VAR))
end

flog=fopen('Selection_parameters.txt','a');
fprintf(flog,'%s\n',datestr(now,'yyyy-mm-ddTHH:MM:SS'));
for k=1:length(NF);
    VAR=getfield(CFG,char(NF(k)));
    fprintf(flog,'%s = %s\n',char(NF(k)),string(VAR));
end
%fclose(flog);
% keyboard

%% Check if NLL DIRECTORY exist
EXI=exist(CFG.NLLOUT);
if EXI==0
    mkdir(CFG.NLLOUT);
end

%keyboard
%% Work
D=dir([CFG.MATDIR '*.mat']);

%% process only non empty files
for j=1:numel(D);
    if D(j).bytes==273 | D(j).bytes==0
        LOG(j)=1;
    else
        LOG(j)=0;
    end
end
i=find(LOG==1);
D(i)=[];

fprintf('Total number of .mat files= %d - non empty %d - removed %d\n',numel(LOG),numel(D),numel(i))


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
                    SLA(L).Punc= S.SLAVE(j).Punc;
                    SLA(L).Pcc=  S.SLAVE(j).Pcc;
                    SLA(L).Ssta= S.SLAVE(j).Ssta;
                    SLA(L).Stim= S.SLAVE(j).Stim;
                    SLA(L).Sunc= S.SLAVE(j).Sunc;
                    SLA(L).Scc=  S.SLAVE(j).Scc;
                    SLA(L).MTP=MTP;
                    SLA(L).MTS=MTS;
                    SLA(L).MP=MP;
                    SLA(L).MS=MS;
                    SLA(L).Pamp=  S.SLAVE(j).Pamp;
                    SLA(L).Pwvf=  S.SLAVE(j).Pwvf;
                    SLA(L).TPamp= S.SLAVE(j).TPamp;
                    SLA(L).TPwvf= S.SLAVE(j).TPwvf;
                    SLA(L).Samp=  S.SLAVE(j).Samp;
                    SLA(L).Swvf=  S.SLAVE(j).Swvf;
                    SLA(L).TSamp= S.SLAVE(j).TSamp;
                    SLA(L).TSwvf= S.SLAVE(j).TSwvf;
                    SLA(L).HEADER=S.SLAVE(j).HEADER;
                    SLA(L).MAGP   =mean(SLA(L).HEADER{6}+log10(SLA(L).Pamp)-log10(SLA(L).TPamp));
                    SLA(L).MAGS   =mean(SLA(L).HEADER{6}+log10(SLA(L).Samp)-log10(SLA(L).TSamp));
                    SLA(L).MAG=mean([SLA(L).MAGP SLA(L).MAGS]);
                    SLA(L).MAGERR=std([SLA(L).MAGP SLA(L).MAGS]);
%                     %% DEBUG
%                     [CFG.MATDIR D(k).name]
%                     SLA(L)
%                     length(SLA(L).Psta)-length(SLA(L).Punc)
%                     length(SLA(L).Ssta)-length(SLA(L).Sunc)
%                     pause
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
fprintf(flog,'Total # of events with CC > %3.2f & NP >= %d & NS >= %d : %d\n',CFG.CCTHR,CFG.NP,CFG.NS,length(SLA));


for k=1:length(SLA)
    PTIME(k)=SLA(k).Ptim(1);
end
[PTIME ID]=sort(PTIME);
SLA=SLA(ID);
T=struct2table(SLA);

% keyboard
%% RIMOZIONE DEI DOPPIONI O DEGLI EVENTI CHE SONO A DISTANZA < CFG.DT
TOL=CFG.DT/86400;
[C,ID,IJ] = uniquetol(PTIME,TOL/max(abs(PTIME)));
NMAX=max(IJ);
TOERASE=[];
%%
for j=1:NMAX
    NSTAP=[];NSTAS=[];CCP=[];CCS=[];Pck=[];
    II=find(IJ==j);
    for k=1:length(II)
        try
        Pck (k)=T{II(k),'Ptim'}{1}(1); %% improvvisamente fa errore
        catch
        Pck (k)=T{II(k),'Ptim'}(1);
        end
        %NSTAP(k)=length(T{II(k),'Psta'}{:});
        %NSTAS(k)=length(T{II(k),'Ssta'}{:});
        NSTAP(k)=length(SLA(II(k)).Psta);
        NSTAS(k)=length(SLA(II(k)).Ssta);
        CCP  (k)=(T{II(k),'MP'});
        CCS  (k)=(T{II(k),'MS'});
    end
    if numel(II)>1  
        % datestr(C(j));
        Z=[NSTAP ; NSTAS; CCP; CCS];
        SCORE=sum(Z);
        [w v]=max(SCORE);
        CHOICE=II(v);       
        OUT=setdiff(II,CHOICE);
        LO=length(OUT);
        LE=length(TOERASE);
        TOERASE(LE+1:LE+length(OUT))=OUT ;
    end
end
SLA(TOERASE)=[];

fprintf('Detected %4.0f double events at first pass (DT <= %3.1f s)\n', numel(TOERASE),CFG.DT)
fprintf(flog,'Detected %4.0f double events at first pass (DT <= %3.1f s)\n', numel(TOERASE),CFG.DT);

%% DOUBLE CHECK 
PTIME=[];    
for k=1:length(SLA)
    PTIME(k)=SLA(k).Ptim(1);
end
[PTIME ID]=sort(PTIME);
SLA=SLA(ID);
T=struct2table(SLA);
TOL=CFG.DT/86400;
[C,ID,IJ] = uniquetol(PTIME,TOL/max(abs(PTIME)));
NMAX=max(IJ);
TOERASE=[];
for j=1:NMAX
    NSTAP=[];NSTAS=[];CCP=[];CCS=[];Pck=[];
    II=find(IJ==j);
    for k=1:length(II)
        try
        Pck (k)=T{II(k),'Ptim'}{1}(1); %% improvvisamente fa errore
        catch
        Pck (k)=T{II(k),'Ptim'}(1);
        end
%        NSTAP(k)=length(T{II(k),'Psta'}{:});
%        NSTAS(k)=length(T{II(k),'Ssta'}{:});
        NSTAP(k)=length(SLA(II(k)).Psta);
        NSTAS(k)=length(SLA(II(k)).Ssta);

        CCP  (k)=(T{II(k),'MP'});
        CCS  (k)=(T{II(k),'MS'});
    end
    if numel(II)>1         
        % datestr(C(j));
        Z=[NSTAP ; NSTAS; CCP; CCS];
        SCORE=sum(Z);
        [w v]=max(SCORE);
        CHOICE=II(v);       
        OUT=setdiff(II,CHOICE);
        LO=length(OUT);
        LE=length(TOERASE);
        TOERASE(LE+1:LE+length(OUT))=OUT ;
    end
end
SLA(TOERASE)=[];
fprintf('Detected %4.0f double events at second pass (DT <= %3.1f s)\n', numel(TOERASE),CFG.DT)
fprintf(flog,'Detected %4.0f double events at second pass (DT <= %3.1f s)\n', numel(TOERASE),CFG.DT);

%keyboard

%% NOW CHEK HOW MANY TEMPLATES ARE PRESENT AND SET SLA.ISTEMPL FIELD ==1
% FIRST SET ISTEMPL=0 FOR ALL the NEW DETECTED EVENTS:
for k=1:numel(SLA)
    SLA(k).ISTEMPL=0;
end

% NOW LOAD TEMPLATES
if CFG.JSON==0
D=dir('../TEMPLATES/*.xml');
else
D=dir('TEMPLATES/*.json');
end
fprintf('Loaded %3.0f templates\n',numel(D));
%% FIND THE TEMPLATES

% keyboard
nT=0;
for k=1:numel(D);
if CFG.JSON==0
    A=mcorr_io_decodeQML(['../TEMPLATES/' D(k).name]);
else
    try
    A=mcorr_io_decodeJSON(['TEMPLATES/' D(k).name]);
    catch
        keyboard
    end
end
    namex=fieldnames(A);
    PARR=[];
    for j=1:length(namex)-1
        XU=getfield(A,char(namex(j)));
        if XU.P > 0
            try
                PARR(j)=XU.P ;
            catch
                keyboard
                %keyboard
            end
        end
    end
    [PARR i]=sort(PARR);
    namex=namex(i);
    
    for n=1:numel(SLA);
        hh=0; ii=[];
        while isempty(ii)==1
            hh=hh+1;
            try
                ii=find(contains(namex,SLA(n).Psta(hh)));            
            catch
                ii=1;
            end
            try
                TD=abs(SLA(n).Ptim(hh)-A.(char(namex(ii))).P)*86400;
            catch
                % fprintf('** \n')
            end
        end
        if TD < 0.5
            nT=nT+1;
            % keyboard
            SLA(n).HEADER=A.HEADER;
            SLA(n).ISTEMPL=1;
            SLA(n).MAGP  = A.HEADER{6};
            SLA(n).MAGS  = A.HEADER{6};
            SLA(n).MAG   = A.HEADER{6};
            SLA(n).MAGERR= 0;       
            TD=[];
        end
    end
end

if CFG.MAGFILT==1
    %% APPLY A FILTER ACCORDING TO MAG ERROR
    % std is calculated only for new detected events
    n=0; MERR=zeros(1,10000);
    for k=1:numel(SLA)
        %    if SLA(k).ISTEMPL==0
        n=n+1;
        MERR(k)=SLA(k).MAGERR;
        %    end
    end
    MERR(n+1:end)=[];
    i=find(MERR > CFG.MULMAG*std(MERR));
    SLA(i)=[];
    fprintf('After magnitude check %d events were removed (mag std: %5.3f)\n',numel(i),std(MERR));
    fprintf(flog,'After magnitude check %d events were removed (mag std: %5.3f)\n',numel(i),std(MERR));
    
    MERR(i)=[];
    fprintf('After magnitude check new mag std: %5.3f\n\n',std(MERR));
    fprintf(flog,'After magnitude check new mag std: %5.3f\n\n',std(MERR));
else
    fprintf('No magnitude filter applied\n');
    fprintf(flog,'No magnitude filter applied\n');
end
%%
SLAVE=SLA;
fprintf('TOTAL NUMBER OF EVENTS COLLECTED: %4d\n',numel(SLAVE))
fprintf('NUMBER OF TEMPLATES FOUND       : %4d\n',nT)
fprintf('NUMBER OF NEW EVENTS FOUND      : %4d\n\n', length(SLAVE)-nT)

% keyboard
fprintf(flog,'TOTAL NUMBER OF EVENTS COLLECTED: %4d\n',numel(SLAVE));
fprintf(flog,'NUMBER OF TEMPLATES FOUND       : %4d\n',nT);
fprintf(flog,'NUMBER OF NEW EVENTS FOUND      : %4d\n\n', length(SLAVE)-nT);
fclose(flog);


%%
% PLOT SERIES
for k=1:length(SLAVE); 
    OTIME(k)=SLAVE(k).TemplateOTime;
    TIMEPP(k)=SLAVE(k).Ptim(1);
    TIMES(k)=SLAVE(k).Stim(1);
    CCP(k)=SLAVE(k).MP; 
end


% EVtime(1)=datenum(2022,5,3,15,50,49);EVM(1)=3.7;
% EVtime(2)=datenum(2022,5,3,20,14,20);EVM(2)=3.5;
% EVtime(3)=datenum(2022,5,10,3,51,17);EVM(3)=3.5;
% EVtime(4)=datenum(2022,5,12,21,12,3);EVM(4)=3.7;

% keyboard
%%
if CFG.FIG==1
    figure;
    subplot(211)
    scatter(TIMEPP,1:length(TIMEPP),30,CCP,'s','filled'); hold on
    datetick('x','mmdd','keeplimits')
    box on
    grid on
    ylabel('Cumulated #')
    colorbar('East')
    title(sprintf('%s - %s - NP=%d NS=%d CC=%4.2f',datestr(floor(TIMEPP(1)),'yyyymmddHHMM'),datestr(ceil(TIMEPP(end)),'yyyymmddHHMM'),CFG.NP,CFG.NS,CFG.CCTHR))
    
    % hold on
    % scatter(EVtime,ones(1,4)*r,EVM.^4,'pk','MarkerFaceColor','y')
    
%    COL=zeros(length(TIMEPP),3);
    COL=zeros(length(TIMEPP),1);
    for k=1:length(TIMEPP)
        COL(k)=SLAVE(k).ISTEMPL;
    end
    %%
    subplot(212);
    i=find(COL==0);
    scatter(TIMEPP(i),i,10,'r','filled')
    hold on
    i=find(COL==1);
    scatter(TIMEPP(i),i,10,'MarkerFaceColor',[.2 .2 .2],'MarkerEdgeColor','none','MarkerFaceAlpha',.4)
    hold off
    
    %plot(TIMEPP,1:length(TIMEPP),'sk','MarkerFaceColor','r','MarkerSize',4);
    datetick('x','mm/dd','keeplimits');grid on
    xlabel('Time ')
    ylabel('Cumulated #')
    box on
    subplot(212);
    n=0;
    for k=1:length(SLAVE)
        MAG(k)=SLAVE(k).MAG;
        MAGE(k)=SLAVE(k).MAGERR;
        if SLAVE(k).ISTEMPL==1
            n=n+1;
            TMAG(n)=SLAVE(k).MAG;
            TMAGE(n)=SLAVE(k).MAGERR;
        end
    end
%    keyboard
    figure;
    %subplot(211)
    EDGES=(floor(min(MAG)):.1:ceil(max(MAG)));
    [N,EDGES] = histcounts(MAG,EDGES);
   % bar(EDGES(1:end-1),log10(fliplr(cumsum(fliplr(N)))),'r')
    plot(EDGES(1:end-1),log10(fliplr(cumsum(fliplr(N)))),'sk','MarkerSize',10,'MarkerFaceColor','r')

    hold on
    EDGES=(floor(min(TMAG)):.1:ceil(max(TMAG)));
    [N,EDGES] = histcounts(TMAG,EDGES);
    %bar(EDGES(1:end-1),log10(fliplr(cumsum(fliplr(N)))),'k')
    plot(EDGES(1:end-1),log10(fliplr(cumsum(fliplr(N)))),'sk','MarkerFaceColor','k')
    grid 
    xlabel('Estimated Magnitude')
    ylabel('Log_1_0 Number')
    legend('Final Catalog','Initial Catalog')
    % keyboard
    figure; %subplot(212)
    IDST=[];
    IDSNO=[];
    for k=1:length(SLAVE); 
        if SLAVE(k).ISTEMPL==1
            IDST=[IDST k];
        else
            IDSNO=[IDSNO k];
        end
    end
        
    errorbar(MAG,MAGE,'.')
    hold on
    plot(MAG,'.k','LineWidth',2);
    % keyboard
    hold on
    plot(IDSNO,MAG(IDSNO),'.','Color',[.9 .0 .0],'LineWidth',2);   
    xlabel('Event ID'); ylabel('Estimated Magnitude'); grid on
    % hold on
    % scatter(EVtime,ones(1,4)*100,EVM.^4,'pk','MarkerFaceColor','y')
    
    figure; %% some other statistical things
    for k=1:length(SLAVE);
        PCC(k)=SLAVE(k).MP;
        SCC(k)=SLAVE(k).MS;
        NP (k)=numel(SLAVE(k).Psta);
        NS (k)=numel(SLAVE(k).Ssta);
        TIM(k)=SLAVE(k).Ptim(1);
        MAG(k)=SLAVE(k).MAG;
    end
    subplot(221); histogram(PCC);                   xlabel('Pcc')     ;ylabel('counts')
    subplot(222); histogram(SCC,'FaceColor','r');   xlabel('Scc')     ;ylabel('counts')
    subplot(223); histogram(NP);                    xlabel('# P phs') ;ylabel('counts')
    subplot(224); histogram(NS,'FaceColor','r');    xlabel('# S phs') ;ylabel('counts')

    
end

% %% BVAL MATTEO BAGAGLI 8/5/2023
% if CFG.FIG==1
%     figure;
%     subplot(111)
% %     [fBValue, fStdDev, fAValue] = bvalue_maxlike2(S.MAG,2.5,0.1, 'DAJE', 'ML');
% %     set(gca,'fontsize',25)
% %     xlim([-2  5])
% %     ylim([10^-0.5 10^5])
%     
%     Mco = min(S.MAG):0.1:max(S.MAG); %min(mag(idx_h9)) cut-off to tests for Mc on a given cluster
%     for i=1:length(Mco)        
%         [fBValue, fStdDev, fAValue] = bvalue_maxlike2( ...
%                                             S.MAG, Mco(i), 0.1, ...                                            'DAJE', 'ML');
%         b_McTest(i)=fBValue;
%         std_McTest(i)=fStdDev;
%     end
% b_valH=    1.0194;  %for deltaM=0.01  %Full1.0194  MQ1.0036
% % figure(4)
% errorbar(Mco(1:length(b_McTest)),b_McTest,std_McTest,'ko','Markersize',8);
% hold on
% plot(Mco(1:length(b_McTest)),ones(1,length(b_McTest)).*b_valH,'--k')
% xlabel('Mco')
% ylabel('b-value' )
% box on
% legend('b-value estimates','correct b-value','Location','Northwest')
% set(gca,'fontsize',25)
% xlim([-1.5 4])  %2.7
% ylim([0 3]) %6.2    
% end

ANS=input('Save RESULTS? [0/1]: ');
if ANS==1
	save('RESULTS.mat','SLAVE','NF')
end
% keyboard
ANS=input('Write NLL OBS? [0/1]: ');

if ANS==1
    %% % WRITE OUT NNLINLOC PHS FILE
    fprintf('Removing old nlloc observation files...\n')
    delete ([CFG.NLLOUT '*'])
    fprintf('Writing new nlloc observation files...\n')
    STHR=sprintf('%3.2f',THR);
    
    for k=1:length(SLAVE)
%         k
%         SLAVE(k)
%         pause
%         SLAVE(k).Psta;
        NAME=[datestr(SLAVE(k).Ptim(1),'yyyymmdd_HHMMSS.FFF') sprintf('_%4.2f',SLAVE(k).MAG)];
        fout=fopen([CFG.NLLOUT NAME '.' STHR '.xml'],'w');
        for j=1:length(SLAVE(k).Psta);
            PSTA=SLAVE(k).Psta(j);
            PCH ='Z';
            PHS ='P';
            TIME=datestr(SLAVE(k).Ptim(j),'yyyymmdd HHMM SS.FFF');
            WT=SLAVE(k).Punc(j);
            WW  =WT+(WT*(1-SLAVE(k).Pcc(j)));
%            WW  =(1-SLAVE(k).Pcc(j))/3;
            %% WT+(WT*(1-CC^2)) \\\--- IMPLEMENTARE QUESTO
            fprintf(fout,'%-6s ?    %s    ? %s      ? %s0 GAU %9.2e -1.00e+00 -1.00e+00 -1.00e+00\n',...
                char(PSTA),PCH,PHS,TIME,WW);
        end
        for j=1:length(SLAVE(k).Ssta);
            SSTA=SLAVE(k).Ssta(j);
            SCH ='Z';
            PHS ='S';
            TIME=datestr(SLAVE(k).Stim(j),'yyyymmdd HHMM SS.FFF');
            WT=SLAVE(k).Sunc(j);
            WW  =WT+(WT*(1-SLAVE(k).Scc(j)));
            %WW  =(1-SLAVE(k).Scc(j))/3;
            fprintf(fout,'%-6s ?    %s    ? %s      ? %s0 GAU %9.2e -1.00e+00 -1.00e+00 -1.00e+00\n',...
                char(SSTA),SCH,PHS,TIME,WW);
        end
        fclose(fout);
    end
end

%%





% %%
% ZIPFILE=sprintf('zip IMPRUNETA_%1dp_%1ds_%4.2f.zip NLL_OBS/*.xml',NP,NS,THR);
% eval(ZIPFILE)

%%
load RESULTS.mat
for NSEIS=1:100;
    ANS=input('Vuoi vedere i sismogrammi? [0/1] ');
    if ANS==1
        fprintf('Loading RESULTS.mat\n')
        N=input('Quale evento vuoi vedere?');
        EV=SLAVE(N)
        np=numel(EV.Psta);
        figure;
        for k=1:np;
            C=.6;
            % subplot(np,1,k)
            plot(k+ normize(EV.Pwvf(k).data).*C,'k'); hold on
            plot(k+normize(EV.TPwvf(k).data).*C,'r');% xlabel('Samples')
            S=sprintf('%s - %3.2f',char(EV.Psta(k)),EV.Pcc(k));
            text(1,k+.3,S)
        end
        title(sprintf('ID: %4d - P-waves - Mag= %4.2f',N,EV.MAG))
        legend('Match','Template')
        ns=numel(EV.Ssta);
        figure;
        for k=1:ns
            % subplot(ns,1,k)
            plot(k+normize(EV.Swvf(k).data).*C,'k'); hold on
            plot(k+normize(EV.TSwvf(k).data).*C,'r'); xlabel('Samples')
            S=sprintf('%s - %3.2f',char(EV.Ssta(k)),EV.Scc(k));
            text(1,k+.3,S)
        end
        title(sprintf('ID_ %4d - S-waves - Mag= %4.2f',N,EV.MAG))
        legend('Match','Template')
    else
        return
    end
end
    
%%    
   



return

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
