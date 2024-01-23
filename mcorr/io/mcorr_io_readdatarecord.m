function [D,offset] = mcorr_io_readdatarecord(f,fid,offset,le,ef,wo,rl,forcebe,verbose,notc,force)
% read_data_record uses global variables f, fid, offset, le, ef, wo, rl,
%   and verbose. It reads a data record and returns a structure D.


fseek(fid,offset,'bof');

% --- read fixed section of Data Header (48 bytes)
D.SequenceNumber        = fread(fid,6,'*char')';
D.DataQualityIndicator  = fread(fid,1,'*char');
D.ReservedByte          = fread(fid,1,'*char');
D.StationIdentifierCode = fread(fid,5,'*char')';
D.LocationIdentifier    = fread(fid,2,'*char')';
D.ChannelIdentifier     = fread(fid,3,'*char')';
D.NetworkCode           = fread(fid,2,'*char')';
D.ChannelFullName = sprintf('%s:%s:%s:%s',deblank(D.NetworkCode), ...
    deblank(D.StationIdentifierCode),deblank(D.LocationIdentifier), ...
    deblank(D.ChannelIdentifier));

% Start Time decoding
[D.RecordStartTime,swapflag] = readbtime(fid,forcebe);
D.RecordStartTimeISO = sprintf('%4d-%03d %02d:%02d:%07.4f',D.RecordStartTime);

if swapflag
    if le
        machinefmt = 'ieee-be';
        le = 0;
    else
        machinefmt = 'ieee-le';
        le = 1;
    end
    position = ftell(fid);
    fclose(fid);
    fid = fopen(f,'rb',machinefmt);
    fseek(fid,position,'bof');
    if verbose > 0
        warning('RDMSEED:DataIntegrity', ...
            'Sequence # %s: need to switch file encoding to %s...\n', ...
            D.SequenceNumber,machinefmt);
    end
end

D.NumberSamples         = fread(fid,1,'uint16');

% Sample Rate decoding
SampleRateFactor        = fread(fid,1,'int16');
SampleRateMultiplier    = fread(fid,1,'int16');
if SampleRateFactor > 0
    if SampleRateMultiplier >= 0
        D.SampleRate = SampleRateFactor*SampleRateMultiplier;
    else
        D.SampleRate = -1*SampleRateFactor/SampleRateMultiplier;
    end
else
    if SampleRateMultiplier >= 0
        D.SampleRate = -1*SampleRateMultiplier/SampleRateFactor;
    else
        D.SampleRate = 1/(SampleRateFactor*SampleRateMultiplier);
    end
end

D.ActivityFlags          = fread(fid,1,'uint8');
D.IOFlags                = fread(fid,1,'uint8');
D.DataQualityFlags       = fread(fid,1,'uint8');
D.NumberBlockettesFollow = fread(fid,1,'uint8');
D.TimeCorrection         = fread(fid,1,'int32');    % Time correction in 0.0001 s
D.OffsetBeginData        = fread(fid,1,'uint16');
D.OffsetFirstBlockette   = fread(fid,1,'uint16');

% --- read the blockettes
OffsetNextBlockette = D.OffsetFirstBlockette;

D.BLOCKETTES = [];
b2000 = 0;  % Number of Blockette 2000

for i = 1:D.NumberBlockettesFollow
    fseek(fid,offset + OffsetNextBlockette,'bof');
    BlocketteType = fread(fid,1,'uint16');

    switch BlocketteType

        case 1000
            % BLOCKETTE 1000 = Data Only SEED (8 bytes)
            OffsetNextBlockette = fread(fid,1,'uint16');
            D.BLOCKETTES.B1000.EncodingFormat = fread(fid,1,'uint8');
            D.BLOCKETTES.B1000.WordOrder = fread(fid,1,'uint8');
            D.BLOCKETTES.B1000.DataRecordLength = fread(fid,1,'uint8');
            D.BLOCKETTES.B1000.Reserved = fread(fid,1,'uint8');

        case 1001
            % BLOCKETTE 1001 = Data Extension (8 bytes)
            OffsetNextBlockette = fread(fid,1,'uint16');
            D.BLOCKETTES.B1001.TimingQuality = fread(fid,1,'uint8');
            D.BLOCKETTES.B1001.Micro_sec = fread(fid,1,'int8');
            D.BLOCKETTES.B1001.Reserved = fread(fid,1,'uint8');
            D.BLOCKETTES.B1001.FrameCount = fread(fid,1,'uint8');

        case 100
            % BLOCKETTE 100 = Sample Rate (12 bytes)
            OffsetNextBlockette = fread(fid,1,'uint16');
            D.BLOCKETTES.B100.ActualSampleRate = fread(fid,1,'float32');
            D.BLOCKETTES.B100.Flags = fread(fid,1,'uint8');
            D.BLOCKETTES.B100.Reserved = fread(fid,1,'uint8');

        case 500
            % BLOCKETTE 500 = Timing (200 bytes)
            OffsetNextBlockette = fread(fid,1,'uint16');
            D.BLOCKETTES.B500.VCOCorrection = fread(fid,1,'float32');
            D.BLOCKETTES.B500.TimeOfException = readbtime(fid,forcebe);
            D.BLOCKETTES.B500.MicroSec = fread(fid,1,'int8');
            D.BLOCKETTES.B500.ReceptionQuality = fread(fid,1,'uint8');
            D.BLOCKETTES.B500.ExceptionCount = fread(fid,1,'uint16');
            D.BLOCKETTES.B500.ExceptionType = fread(fid,16,'*char')';
            D.BLOCKETTES.B500.ClockModel = fread(fid,32,'*char')';
            D.BLOCKETTES.B500.ClockStatus = fread(fid,128,'*char')';

        case 2000
            % BLOCKETTE 2000 = Opaque Data (variable length)
            b2000 = b2000 + 1;
            OffsetNextBlockette = fread(fid,1,'uint16');
            BlocketteLength = fread(fid,1,'uint16');
            OffsetOpaqueData = fread(fid,1,'uint16');
            D.BLOCKETTES.B2000(b2000).RecordNumber = fread(fid,1,'uint32');
            D.BLOCKETTES.B2000(b2000).DataWordOrder = fread(fid,1,'uint8');
            D.BLOCKETTES.B2000(b2000).Flags = fread(fid,1,'uint8');
            NumberHeaderFields = fread(fid,1,'uint8');
            HeaderFields = splitfield(fread(fid,OffsetOpaqueData-15,'*char')','~');
            D.BLOCKETTES.B2000(b2000).HeaderFields = HeaderFields(1:NumberHeaderFields);
            % Opaque data are stored as a single char string, but must be
            % decoded using appropriate format (e.g., Quanterra Q330)
            D.BLOCKETTES.B2000(b2000).OpaqueData = fread(fid,BlocketteLength-OffsetOpaqueData,'*char')';

        otherwise
            OffsetNextBlockette = fread(fid,1,'uint16');

            if verbose > 0
                warning('RDMSEED:UnknownBlockette', ...
                    'Unknown Blockette number %d (%s)!\n', ...
                    BlocketteType,D.ChannelFullName);
            end
    end
end

% --- read the data stream
fseek(fid,offset + D.OffsetBeginData,'bof');

if ~force && isfield(D.BLOCKETTES,'B1000')
    EncodingFormat = D.BLOCKETTES.B1000.EncodingFormat;
    WordOrder = D.BLOCKETTES.B1000.WordOrder;
    D.DataRecordSize = 2^D.BLOCKETTES.B1000.DataRecordLength;
else
    EncodingFormat = ef;
    WordOrder = wo;
    D.DataRecordSize = rl;
end

uncoded = 0;

D.d = NaN;
D.t = NaN;

switch EncodingFormat

    case 0
        % --- decoding format: ASCII text
        D.EncodingFormatName = {'ASCII'};
        D.d = fread(fid,D.DataRecordSize - D.OffsetBeginData,'*char')';

    case 1
        % --- decoding format: 16-bit integers
        D.EncodingFormatName = {'INT16'};
        dd = fread(fid,ceil((D.DataRecordSize - D.OffsetBeginData)/2),'*int16');
        if xor(~WordOrder,le)
            dd = swapbytes(dd);
        end
        D.d = dd(1:D.NumberSamples);

    case 2
        % --- decoding format: 24-bit integers
        D.EncodingFormatName = {'INT24'};
        dd = fread(fid,ceil((D.DataRecordSize - D.OffsetBeginData)/3),'bit24=>int32');
        if xor(~WordOrder,le)
            dd = swapbytes(dd);
        end
        D.d = dd(1:D.NumberSamples);

    case 3
        % --- decoding format: 32-bit integers
        D.EncodingFormatName = {'INT32'};
        dd = fread(fid,ceil((D.DataRecordSize - D.OffsetBeginData)/4),'*int32');
        if xor(~WordOrder,le)
            dd = swapbytes(dd);
        end
        D.d = dd(1:D.NumberSamples);

    case 4
        % --- decoding format: IEEE floating point
        D.EncodingFormatName = {'FLOAT32'};
        dd = fread(fid,ceil((D.DataRecordSize - D.OffsetBeginData)/4),'*float');
        if xor(~WordOrder,le)
            dd = swapbytes(dd);
        end
        D.d = dd(1:D.NumberSamples);

    case 5
        % --- decoding format: IEEE double precision floating point
        D.EncodingFormatName = {'FLOAT64'};
        dd = fread(fid,ceil((D.DataRecordSize - D.OffsetBeginData)/8),'*double');
        if xor(~WordOrder,le)
            dd = swapbytes(dd);
        end
        D.d = dd(1:D.NumberSamples);

    case {10,11,19}
        % --- decoding formats: STEIM-1 and STEIM-2 compression
        %    (c) Joseph M. Steim, Quanterra Inc., 1994
        steim = find(EncodingFormat==[10,11,19]);
        D.EncodingFormatName = {sprintf('STEIM%d',steim)};

        % Steim compression decoding strategy optimized for Matlab
        % -- by F. Beauducel, October 2010 --
        %
        %   1. loads all data into a single 16xM uint32 array
        %   2. gets all nibbles from the first row splitted into 2-bit values
        %   3. for each possible nibble value, selects (find) and decodes
        %      (bitsplit) all the corresponding words, and stores results
        %      in a 4xN (STEIM1) or 7xN (STEIM2) array previously filled with
        %      NaN's. For STEIM2 with nibbles 2 or 3, decodes also dnib values
        %      (first 2-bit of the word)
        %   5. reduces this array with non-NaN values only
        %   6. integrates with cumsum
        %
        % This method is about 30 times faster than a 'C-like' loops coding...

        frame32 = fread(fid,[16,(D.DataRecordSize - D.OffsetBeginData)/64],'*uint32');
        if xor(~WordOrder,le)
            frame32 = swapbytes(frame32);
        end

        % specific processes for STEIM-3
        if steim == 3
            % first bit = 1 means second differences
            SecondDiff = bitshift(frame32(1,:),-31);
            % checks for "squeezed flag"... and replaces frame32(1,:)
            squeezed = bitand(bitshift(frame32(1,:),-24),127);
            k = find(bitget(squeezed,7));
            if ~isempty(k)
                moredata24 = bitand(frame32(1,k),16777215);
                k = find(squeezed == 80);   % upper nibble 8-bit = 0x50
                if ~isempty(k)
                    frame32(1,k) = hex2dec('15555555');
                end
                k = find(squeezed == 96);   % upper nibble 8-bit = 0x60
                if ~isempty(k)
                    frame32(1,k) = hex2dec('2aaaaaaa');
                end
                k = find(squeezed == 112);  % upper nibble 8-bit = 0x70
                if ~isempty(k)
                    frame32(1,k) = hex2dec('3fffffff');
                end
            end
        end

        % nibbles is an array of the same size as frame32...
        nibbles = bitand(bitshift(repmat(frame32(1,:),16,1),repmat(-30:2:0,size(frame32,2),1)'),3);
        x0 = bitsign(frame32(2,1),32);  % forward integration constant
        xn = bitsign(frame32(3,1),32);  % reverse integration constant

        switch steim

        case 1
            % STEIM-1: 3 cases following the nibbles
            ddd = NaN*ones(4,numel(frame32));   % initiates array with NaN
            k = find(nibbles == 1);         % nibble = 1 : four 8-bit differences
            if ~isempty(k)
                ddd(1:4,k) = bitsplit(frame32(k),32,8);
            end
            k = find(nibbles == 2);         % nibble = 2 : two 16-bit differences
            if ~isempty(k)
                ddd(1:2,k) = bitsplit(frame32(k),32,16);
            end
            k = find(nibbles == 3);         % nibble = 3 : one 32-bit difference
            if ~isempty(k)
                ddd(1,k) = bitsign(frame32(k),32);
            end

        case 2
            % STEIM-2: 7 cases following the nibbles and dnib
            ddd = NaN*ones(7,numel(frame32));   % initiates array with NaN
            k = find(nibbles == 1);         % nibble = 1 : four 8-bit differences
            if ~isempty(k)
                ddd(1:4,k) = bitsplit(frame32(k),32,8);
            end
            k = find(nibbles == 2);         % nibble = 2 : must look in dnib
            if ~isempty(k)
                dnib = bitshift(frame32(k),-30);
                kk = k(dnib == 1);      % dnib = 1 : one 30-bit difference
                if ~isempty(kk)
                    ddd(1,kk) = bitsign(frame32(kk),30);
                end
                kk = k(dnib == 2);      % dnib = 2 : two 15-bit differences
                if ~isempty(kk)
                    ddd(1:2,kk) = bitsplit(frame32(kk),30,15);
                end
                kk = k(dnib == 3);      % dnib = 3 : three 10-bit differences
                if ~isempty(kk)
                    ddd(1:3,kk) = bitsplit(frame32(kk),30,10);
                end
            end
            k = find(nibbles == 3);             % nibble = 3 : must look in dnib
            if ~isempty(k)
                dnib = bitshift(frame32(k),-30);
                kk = k(dnib == 0);      % dnib = 0 : five 6-bit difference
                if ~isempty(kk)
                    ddd(1:5,kk) = bitsplit(frame32(kk),30,6);
                end
                kk = k(dnib == 1);      % dnib = 1 : six 5-bit differences
                if ~isempty(kk)
                    ddd(1:6,kk) = bitsplit(frame32(kk),30,5);
                end
                kk = k(dnib == 2);      % dnib = 2 : seven 4-bit differences (28 bits!)
                if ~isempty(kk)
                    ddd(1:7,kk) = bitsplit(frame32(kk),28,4);
                end
            end

        case 3  % *** STEIM-3 DECODING IS ALPHA AND UNTESTED ***
            % STEIM-3: 7 cases following the nibbles
            ddd = NaN*ones(9,numel(frame32));   % initiates array with NaN
            k = find(nibbles == 0);             % nibble = 0 : two 16-bit differences
            if ~isempty(k)
                ddd(1:2,k) = bitsplit(frame32(k),32,16);
            end
            k = find(nibbles == 1);             % nibble = 1 : four 8-bit differences
            if ~isempty(k)
                ddd(1:4,k) = bitsplit(frame32(k),32,8);
            end
            k = find(nibbles == 2);             % nibble = 2 : must look even dnib
            if ~isempty(k)
                dnib2 = bitshift(frame32(k(2:2:end)),-30);
                w60 = bitand(frame32(k(2:2:end)),1073741823) ...
                    + bitshift(bitand(frame32(k(1:2:end)),1073741823),30);  % concatenates two 30-bit words
                kk = find(dnib2 == 0);      % dnib = 0: five 12-bit differences (60 bits)
                if ~isempty(kk)
                    ddd(1:5,k(2*kk)) = bitsplit(w60,60,12);
                end
                kk = find(dnib2 == 1);      % dnib = 1: three 20-bit differences (60 bits)
                if ~isempty(kk)
                    ddd(1:3,k(2*kk)) = bitsplit(w60,60,20);
                end
            end
            k = find(nibbles == 3);             % nibble = 3 : must look 3rd bit
            if ~isempty(k)
                dnib = bitshift(frame32(k),-27);
                kk = k(dnib == 24);     % dnib = 11000 : nine 3-bit differences (27 bits)
                if ~isempty(kk)
                    ddd(1:9,kk) = bitsplit(frame32(kk),27,3);
                end
                kk = k(dnib == 25);     % dnib = 11001 : Not A Difference
                if ~isempty(kk)
                    ddd(1,kk) = bitsign(frame32(kk),27);
                end
                kk = k(dnib > 27);      % dnib = 111.. : 29-bit sample (29 bits)
                if ~isempty(kk)
                    ddd(1,kk) = bitsign(frame32(kk),29);
                end
            end
        end

        % Little-endian coding: needs to swap bytes
        if ~WordOrder
            ddd = flipud(ddd);
        end
        dd = ddd(~isnan(ddd));      % reduces initial array ddd: dd is non-NaN values of ddd

        % controls the number of samples
        if numel(dd) ~= D.NumberSamples
            if verbose > 1
                warning('RDMSEED:DataIntegrity','Problem in %s sequence # %s [%d-%03d %02d:%02d:%07.4f]: number of samples in header (%d) does not equal data (%d).\n', ...
                    D.EncodingFormatName{:},D.SequenceNumber,D.RecordStartTimeISO,D.NumberSamples,numel(dd));
            end
            if numel(dd) < D.NumberSamples
                D.NumberSamples = numel(dd);
            end
        end

        % rebuilds the data vector by integrating the differences
        D.d = cumsum([x0;dd(2:D.NumberSamples)]);

        % controls data integrity...
        if D.d(end) ~= xn
            warning('RDMSEED:DataIntegrity','Problem in %s sequence # %s [%s]: data integrity check failed, last_data=%d, Xn=%d.\n', ...
                D.EncodingFormatName{:},D.SequenceNumber,D.RecordStartTimeISO,D.d(end),xn);
        end

        if D.NumberSamples == 0
            D.d = nan(0,1);
        end

        % for debug purpose...
        if verbose > 2
            D.dd = dd;
            D.nibbles = nibbles;
            D.x0 = x0;
            D.xn = xn;
        end

    case 12
        % --- decoding format: GEOSCOPE multiplexed 24-bit integer
        D.EncodingFormatName = {'GEOSCOPE24'};
        dd = fread(fid,(D.DataRecordSize - D.OffsetBeginData)/3,'bit24=>double');
        if xor(~WordOrder,le)
            dd = swapbytes(dd);
        end
        D.d = dd(1:D.NumberSamples);

    case {13,14}
        % --- decoding format: GEOSCOPE multiplexed 16/3 and 16/4 bit gain ranged
        %   (13): 16/3-bit (bit 15 is unused)
        %   (14): 16/4-bit
        %   bits 15-12 = 3 or 4-bit gain exponent (positive)
        %   bits 11-0 = 12-bit mantissa (positive)
        %   => data = (mantissa - 2048) / 2^gain
        geoscope = 7 + 8*(EncodingFormat==14); % mask for gain exponent
        D.EncodingFormatName = {sprintf('GEOSCOPE16-%d',EncodingFormat-10)};
        dd = fread(fid,(D.DataRecordSize - D.OffsetBeginData)/2,'*uint16');
        if xor(~WordOrder,le)
            dd = swapbytes(dd);
        end
        dd = (double(bitand(dd,2^12-1))-2^11)./2.^double(bitand(bitshift(dd,-12),geoscope));
        D.d = dd(1:D.NumberSamples);

    case 15
        % --- decoding format: US National Network compression
        D.EncodingFormatName = {'USNN'};
        uncoded = 1;

    case 16
        % --- decoding format: CDSN 16-bit gain ranged
        D.EncodingFormatName = {'CDSN'};
        uncoded = 1;

    case 17
        % --- decoding format: Graefenberg 16-bit gain ranged
        D.EncodingFormatName = {'GRAEFENBERG'};
        uncoded = 1;

    case 18
        % --- decoding format: IPG - Strasbourg 16-bit gain ranged
        D.EncodingFormatName = {'IPGS'};
        uncoded = 1;

    case 30
        % --- decoding format: SRO format
        D.EncodingFormatName = {'SRO'};
        uncoded = 1;

    case 31
        % --- decoding format: HGLP format
        D.EncodingFormatName = {'HGLP'};
        uncoded = 1;

    case 32
        % --- decoding format: DWWSSN gain ranged format
        D.EncodingFormatName = {'DWWSSN'};
        uncoded = 1;

    case 33
        % --- decoding format: RSTN 16-bit gain ranged
        D.EncodingFormatName = {'RSTN'};
        uncoded = 1;

    otherwise
        D.EncodingFormatName = {sprintf('** Unknown (%d) **',EncodingFormat)};
        uncoded = 1;

end

if uncoded
    error('Sorry, the encoding format "%s" is not yet implemented.',D.EncodingFormatName);
end

% Applies time correction (if needed)
D.RecordStartTimeMATLAB = datenum(double([D.RecordStartTime(1),0,D.RecordStartTime(2:5)])) ...
    + (~notc & bitand(D.ActivityFlags,2) == 0)*D.TimeCorrection/1e4/86400;
tv = datevec(D.RecordStartTimeMATLAB);
doy = datenum(tv(1:3)) - datenum(tv(1),1,0);
D.RecordStartTime = [tv(1),doy,tv(4:5),round(tv(6)*1e4)/1e4];
D.RecordStartTimeISO = sprintf('%4d-%03d %02d:%02d:%07.4f',D.RecordStartTime);

D.t = D.RecordStartTimeMATLAB;

% makes the time vector and applies time correction (if needed)
if EncodingFormat > 0
    D.t = D.t + (0:(D.NumberSamples-1))'/(D.SampleRate*86400);
end


offset = ftell(fid);
fread(fid,1,'char');    % this is to force EOF=1 on last record.
if feof(fid)
    offset = -1;
end
