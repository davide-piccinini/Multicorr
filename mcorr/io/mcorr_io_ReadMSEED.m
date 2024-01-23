function OUT = mcorr_io_ReadMSEED(filein)

% Read mseed files
% OUT=ReadMSEED(filein)
% OUT is a structure

        %keyboard
        %fprintf('Reading %s\n',filein)
        AA=rdmseed(filein);
        OUT.network=AA(1).NetworkCode;
        OUT.station=AA(1).StationIdentifierCode;
        OUT.channel=AA(1).ChannelIdentifier;
        OUT.dataquality='';
        OUT.location='';
        OUT.type='';
        OUT.sampleRate=AA(1).SampleRate;
        OUT.sampleType='i';
        DIFF=AA(1).RecordStartTimeMATLAB - round(AA(1).RecordStartTimeMATLAB);
        DIFF=abs(DIFF/86400);
        if abs(DIFF) < 1/OUT.sampleRate/86400;
            AA(1).RecordStartTimeMATLAB=round(AA(1).RecordStartTimeMATLAB);
        end
        OUT.dateTimeString=datestr(AA(1).RecordStartTimeMATLAB,'yyyy/mm/dd HH:MM:SS.FFF');
        OUT.dateTime.year=str2num(datestr(AA(1).RecordStartTimeMATLAB,'yyyy'));
        OUT.dateTime.year=str2num(datestr(AA(1).RecordStartTimeMATLAB,'mm'));
        OUT.dateTime.year=str2num(datestr(AA(1).RecordStartTimeMATLAB,'dd'));
        OUT.dateTime.year=str2num(datestr(AA(1).RecordStartTimeMATLAB,'HH'));
        OUT.dateTime.year=str2num(datestr(AA(1).RecordStartTimeMATLAB,'MM'));
        OUT.dateTime.year=str2num(datestr(AA(1).RecordStartTimeMATLAB,'SS.FFF'));
        DL=0;
        for k=1:length(AA)
            DL=DL+length(AA(k).d);
            LEN(k)=AA(k).NumberSamples;
        end
        START=AA(1).t(1);
        END  =AA(end).t(end);
        OUT.sampleCount=DL;
        OUT.numberOfSample=DL;
        OUT.matlabTimeVector=linspace(START,END,DL)';
        OUT.data=vertcat(AA.d);
