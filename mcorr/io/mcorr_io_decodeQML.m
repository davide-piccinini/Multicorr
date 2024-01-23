function X = mcorr_io_decodeQML(filexml,PRT)

% usage : X=decodeQML(filexml,PRT)
% read QuakeML files

% revision dec 17 2021
% revision may 11 2022 : add 5 letter stations case
% revisione may 14 corrected LAT and LON reading
% rev june 22 corrected decodeQML to read NLLoc
% rev oct 25 2022 corrected multiple magnitude issue
% Davide Piccinini

if nargin==1
    PRT=0;
end
%%
fid=fopen(filexml);
CC=textscan(fid,'%s','delimiter','\n');
fclose(fid);
KK=CC{1};

X=[];
n=0;

%% EventID
%% PREFERRED ORIGIN
EVID=find(contains(KK,'eventId'));
i=findstr(char(KK(EVID)),'=');
f=findstr(char(KK(EVID)),'">');
L=char(KK(EVID));
EVIDs=L(i(2)+1:f-1);
EVID=str2num(EVIDs);

%% PREFERRED ORIGIN
PREFORIGIN=find(contains(KK,'<preferredOriginID>'));
i=findstr(char(KK(PREFORIGIN)),'=');
f=findstr(char(KK(PREFORIGIN)),'</pref');
L=char(KK(PREFORIGIN));
PREFORIGID=L(i+1:f-1);
i=find(contains(KK,PREFORIGID));
i=i(2);

OTime=char(KK(i+2)); OTime=OTime(8:end-12);MOTime=datenum(OTime,'yyyy-mm-ddTHH:MM:SS.FFF');
OLAT=char(KK(i+6));fi=findstr(OLAT,'</value>');LAT=OLAT(8:fi-1);LAT=str2double(LAT);
OLON=char(KK(i+10));fi=findstr(OLON,'</value>');LON=OLON(8:fi-1);LON=str2double(LON);
ODEP=char(KK(i+14)); fi=findstr(ODEP,'</value>');DEP=ODEP(8:fi-1);DEP=str2double(DEP)/1000;

%% PREFERRED MAGNITUDE
PREFMAG=find(contains(KK,'<preferredMagnitudeID>'));
i=findstr(char(KK(PREFMAG)),'=');
f=findstr(char(KK(PREFMAG)),'</pref');
L=char(KK(PREFMAG));
PREFMAGID=L(i+1:f-1);
i=find(contains(KK,PREFMAGID));
i=i(2);
OMAG=char(KK(i+2));MAG=OMAG(8:10);MAG=str2double(MAG);


% % DECODE LOCATION
% for i=1:length(KK)
%     LIN=char(KK(i));
%     if findstr(LIN,'<origin publicID')==1
%         OTime=char(KK(i+2));
%         OTime=OTime(8:end-12);
%         MOTime=datenum(OTime,'yyyy-mm-ddTHH:MM:SS.FFF');
%     end
%     if findstr(LIN,'<latitude')==1
%         OLAT=char(KK(i+1));fi=findstr(OLAT,'</value>');
%         LAT=OLAT(8:fi-1);LAT=str2double(LAT);
%     end
%     if findstr(LIN,'<longitude')==1
%         OLON=char(KK(i+1));fi=findstr(OLON,'</value>');
%         LON=OLON(8:fi-1);LON=str2double(LON);
%     end
%     if findstr(LIN,'<depth>')==1
%         ODEP=char(KK(i+1)); fi=findstr(ODEP,'</value>');
%         DEP=ODEP(8:fi-1);DEP=str2double(DEP)/1000;
%     end
%     if findstr(LIN,'<author>Sala Sismica INGV-Roma</author>')==1 %%% to account only for revised magnitude
%         OMAG=char(KK(i-9));MAG=OMAG(8:10);MAG=str2double(MAG);
%     end
% end

if exist('MAG')==0; MAG=9; end

%%
% DECODE ARRIVALS
% n=0;
% for i=1:length(KK)
%     LIN=char(KK(i));
%     if findstr(LIN,'arrival publicID')==2
%         n=n+1;
%         PickID(n)=(i+1);
%     end
% end
% fprintf('Number of arrivals: %d\n',n);
%
%%
INI=find(contains(KK,'arrival publicID'));
FIN=find(contains(KK,'</arrival>'));
if PRT==1
fprintf('Number of arrivals: %d\n',length(INI));
end

%%
for i=1:length(FIN)
    BLK=KK(INI(i):FIN(i));
    j=find(contains(BLK,'<timeWeight>'));
    CC=char(BLK(j));
    w =[findstr(CC,'>') findstr(CC,'<')];
    TW(i)=str2double(CC(w(1)+1:w(4)-1));
end
%%

IND=find(TW > 0);
if PRT==1
fprintf('Number of valid picks: %d\n',length(IND));
end

for k=1:length(IND)
    it=IND(k);
    BLK=KK(INI(it):FIN(it));
    i=find(contains(BLK,'pickID'));
    LL=char(BLK(i));
    f=findstr(LL,'=');
    ID=LL(f+1:end-9);
    PICKID{k}=ID;
    if isempty(PICKID{k})==1
        f=findstr(LL,':');
        ID=LL(f+7:end-9);
        PICKID{k}=ID;
    end
end

% fprintf(' Total number of picks: %d\n',length(PICKID))
%%
NENDPICK=find(contains(KK,'</pick>'));
    for k=1:length(PICKID)
        PID=PICKID{k};
        NINI=find(contains(KK,PID));
        BLKI(k)=NINI(2);
        DJ=NENDPICK-BLKI(k);
        i=find(DJ > 0);J=(i(1));
        BLKF(k)=NENDPICK(J);
     end


%%
sn=0;pn=0;
for k=1:length(PICKID)
    PICKID{k};
    BLK=KK(BLKI(k):BLKF(k));
    %pause
    i=find(contains(BLK,'stationCode'));
    LL=char(BLK(i));
    ast=findstr(LL,'"');
    staz=LL(ast(end-1):ast(end)); staz=staz(2:end-1);
    if isfield(X,staz)==0
       X.(staz)  =[];
       X.(staz).P=[];
       X.(staz).Punc=[];
       X.(staz).S=[];
       X.(staz).Sunc=[];
    end

    LL=char(BLK(3));TIMEPICK=LL(8:32) ; TIME    =datenum(TIMEPICK,'yyyy-mm-ddTHH:MM:SS.FFF');
    LL=char(BLK(4));UNC     =LL(14:16); UNC=str2double(UNC);

    i=find(contains(BLK,'phaseHint'));
    LL=char(BLK(i));
    PHASE=LL(12);
    if strcmp(PHASE,'P')
        X.(staz).P=TIME;
        X.(staz).Punc=UNC;
        pn=pn+1;
    end
    if strcmp(PHASE,'S')
        X.(staz).S=TIME;
        X.(staz).Sunc=UNC;
        sn=sn+1;
    end

end


%%
if isempty(X)==0
    X.HEADER={OTime,MOTime,LAT,LON,DEP,MAG,EVID};
end

if PRT==1
    NM=fieldnames(X);
    NTOT=length(NM);
    for k=1:NTOT-1
        PT=X.(char(NM(k))).P;
        ST=X.(char(NM(k))).S;
        fprintf('%5s %s %s \n', char(NM(k)),datestr(PT,'yyyy-mm-ddTHH:MM:SS.FFF'),datestr(ST,'HH:MM:SS.FFF'))
    end
end
