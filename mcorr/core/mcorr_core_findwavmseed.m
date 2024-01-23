function []=mcorr_core_findwavmseed(fileday,filesamp,THR,PHASE);

%% usage: []=findwav_mseed(fileday,filesamp,THR,PHASE);
%% Based on Matching filter algorithm
%
%  filesamp = file sac che contiene la f.o. campione
%  fileday  = file sac nel quale bisogna effettuare la ricerca
%  fvar     = 1 se il filesamp ? una variabile o un file fisico [optional]
%
% Controllare la durata del match (DUR)
% la larghezza del filtro
% e il valore della CC
%
% Crea in uscita un file "wavl.out"
%

% implementato riducendo considerevolmente
% la finestra di analisi

% Aggiunto il dimensionamento iniziale della variabile MC

% Mod. 16/07/08
% riordinata la parte della scrittura
% commentato alcune i blocchi principali
% aggiunto il sort delle var per i files di uscita

% Mod 31/10/2020
% i files sono entrambi in memoria


PRT=0;

DIV=4;

DAY=fileday;
DAYNAME=datestr(datenum(DAY.dateTimeString,'yyyy/mm/dd HH:MM:SS.FFF'),'yyyymmdd');

NAMEOUT=[filesamp.OTime '.' filesamp.station '.' filesamp.channel '.' PHASE '.' DAYNAME '.out'];

if strcmpi(PHASE,'P')==1
    LW=length(filesamp.ptemplate);
    if (LW/DIV)-floor(LW/DIV) > 0   %% SE DISPARI
        WAVLET=filesamp.ptemplate;
        WAVLET=WAVLET(1:end-1);
        LW=length(WAVLET);
    else
        WAVLET=filesamp.ptemplate;
        LW=length(WAVLET);
    end
else
    LW=length(filesamp.stemplate);
    if (LW/DIV)-floor(LW/DIV) > 0   %% SE DISPARI
        WAVLET=filesamp.stemplate;
        WAVLET=WAVLET(1:end-1);
        LW=length(WAVLET);
    else
        WAVLET=filesamp.stemplate;
        LW=length(WAVLET);
    end
end

LD=length(fileday.datafilt);
SHIFT=floor(LW/DIV);
NWIN=floor((LD-LW)/SHIFT)+1;
%% inizializza CCMAX !! <<-- questo fa guardagnare un sacco di tempo
CCMAX=single(zeros(1,NWIN));
%% INIZIALIZZA OUTSAC e STARTORI
STARTORI=single(zeros(1,NWIN));
L=length(DAY.datafilt);
NW=ceil(L/LW);

% disp(sprintf('*** CC THRESHOLD = %3.2f',THR))
if PRT==1
fprintf('Matching started @ %s\n',DAY.station)
end
t0 = clock;
N=0;

TIM=single(zeros(1,NWIN));
MAXTRAX=zeros(1,NWIN);
% keyboard
%% original code
for J=1:NWIN;
    INI=(J-1)*SHIFT+1;
    FIN=INI+LW-1;
    TRAX=DAY.datafilt(INI:FIN);
    [C,~]=xcorr(WAVLET,TRAX,'coeff');
    %[C,~]=xcorr(zscore(WAVLET),zscore(TRAX),'unbiased');
    [CMAX,JMAX]=max((C));  %% ABS o non ABS?  Non ABS perché deve cercare SOLO queli identici e non con polarità invertita
    if CMAX >= THR %% supera la soglia della CC --> FILTER MATCHED
        TMAX=round(LW-JMAX);
        I1=INI+TMAX; if I1 < 1; I1=1; end
        I2=I1 +LW-1;
        TIM(J)=I1;      % PRIMO CAMPIONE DEL MATCH
        if TIM(J) > 0
            try
                MAXTRAX(J)=max(DAY.datafilt(I1:I2)); %% allineato sullo stack
            catch
                fprintf('%s %s -->> %d %d %d --- %d / %d \n',filesamp.station,PHASE,I1,I2,length(DAY.datafilt),J,NWIN)
                MAXTRAX(J)=max(DAY.datafilt(I1:length(DAY.datafilt)));
            end
            CCMAX  (J)=CMAX;
            N=N+1;
            STARTORI(N)=I1-1;
        end
    end
end
% disp(sprintf('...done in %d secs!',round(etime(clock,t0))));
if PRT==1
fprintf('Corr @ %s for %s done in %5.0f secs!\n',DAY.station,PHASE,round(etime(clock,t0)));
end
%disp(' ');
%return
%%
if isempty (MAXTRAX)==0
    %% ESTRAE SOLO I DATI  DHE HANNO MATCHATO IL FILTRO %%            
    IDAMP=find(MAXTRAX > 0);% indice in base alle ampiezze
    CCOUT=CCMAX(IDAMP);
    ITIM =TIM(IDAMP);
%     CCOUT=CCMAX(find(CCMAX > 0));  %MAX della CC
%     ITIM =TIM(find(TIM > 0));
    IDAMP=IDAMP';
    MT   =MAXTRAX(IDAMP)';
    %keyboard
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PER EVITARE CHE LO STESSO EVENTO SIA STATO MATCHATO 2 VOLTE
    %% CONSECUTIVAMENTE, CONTROLLO LE AMPIEZZE E TOLGO LE COPPIE DI VALORI
    %% UGUALI
    [~,I,J]=unique(MT,'legacy');
    %% ROUTINE PER SCEGLIERE I IN FUNZIONE DEL VALORE DI CC PIU' ALTO
    for k=1:max(J);
        i=find(J==k);
        if length(i)>1
            [~, w]=max(CCOUT(i));
            igood=i(w);
            I(k)=igood;
        end
    end
    %%
    NMT=MT(I);
    NCCOUT=CCOUT(I);
    NTIM   =ITIM(I);
    NSTARTORI=STARTORI(I);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SORT DELLE USCITE IN FUNZIONE DEL TEMPO
    [~,J] = sort(NTIM);
    SNMT      = NMT(J);
    SNCCOUT   = NCCOUT(J);
    SNSTARTORI= NSTARTORI(J);
    %fprintf('---->>> %d\n',numel(J))
    try
    if numel(J) > 0 && SNSTARTORI(1)==0; SNSTARTORI(1)=1;end %% SE IL PRIMO CAMPIONE E' 0
    catch
    fprintf('Houston we had a problem here in findwav...')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    fout=fopen(['WORK/' NAMEOUT],'w');
    for i=1:length(SNCCOUT);
        OUTTIME=datestr(DAY.matlabTimeVector(SNSTARTORI(i)),'yyyy-mm-ddTHH:MM:SS.FFF');
        fprintf(fout,'%s %8.6e %5.3f\n',OUTTIME,SNMT(i), SNCCOUT(i));%
    end
    fclose (fout);
end

%clear fileday filesamp SNSTARTORI SNCCOUT SNMT NTIM CCMAX STARTORI TIM DAY NMT NCCOUT NTIM NSTARTORI
clear all

return
