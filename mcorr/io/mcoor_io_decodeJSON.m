function X=mcorr_io_decodeJSON(filein)

% X=mcoor_io_decodeJSON(filein)

str=fileread(filein);
data = jsondecode(str);

FMT ='yyyy-mm-ddTHH:MM:SS.FFF';
FMTM='yyyy-mm-dd HH:MM:SS.FFF';
MOTIME=datenum(data.date,'yyyy-mm-dd HH:MM:SS.FFF');
OTIME=datestr(MOTIME,FMT);

ID=data.id;
LAT=data.location.latitude;
LON=data.location.longitude;
DEP=data.location.depth;
MAG=data.magnitudos.value;

X.HEADER={OTIME,MOTIME,LAT,LON,DEP,MAG,ID};

N=numel(fieldnames(data.location.picks,'-full'));
STLIST=fieldnames(data.location.picks,'-full');

for k=1:N
    STA=STLIST{k};
    try
        PTIME=getfield(data.location.picks.(STA).P,'pick_time');
        Punc=0.05;
        X.(STLIST{k}).P=datenum(PTIME,FMTM);
        X.(STLIST{k}).Punc=Punc;
    catch
        %fprintf('NO P for %s \n',STA);
        X.(STLIST{k}).P=[];
        X.(STLIST{k}).Punc=[];
    end
    try
        STIME=getfield(data.location.picks.(STA).S,'pick_time');
        Sunc=0.05;
        X.(STLIST{k}).S=datenum(STIME,FMTM);
        X.(STLIST{k}).Sunc=Sunc;

    catch
        %fprintf('NO S for %s \n',STA);
        X.(STLIST{k}).S=[];
        X.(STLIST{k}).Sunc=[];
    end
    
end
