function FULLDAY=mcorr_io_readDAY(DAYDIR,FILT)

% read in all the files in DAYDIR
% and store in the FULLDAY structure
% could be parallelized

FILT;
FULLDAY={};

D=dir([DAYDIR '/*.*.*HZ*']);
n=numel(D);
fprintf('Reading from %s [%2d]\n',DAYDIR,n)

parfor i=1:n
%for i=1:numel(D)    
    DAY=mcorr_io_ReadMSEED([D(i).folder '/' D(i).name]);
    trace=DAY.data;
    DAY.DAYDIR=DAYDIR;
    if length(trace) > 1000  %% skip dayfiles with no data
       trace=trace-mean(trace);
       trace=mcorr_signal_taper(trace,100);
       DAY.datafilt = mcorr_signal_filtrax(trace,FILT(1),FILT(2),DAY.sampleRate);
       DAY.DAYDIR=DAYDIR;
    end
    FULLDAY{i}=DAY;
end
S='';
for k=1:n
    S=[S sprintf('%s.%s',FULLDAY{k}.network,FULLDAY{k}.station)];
end
fprintf('%s\n',S);