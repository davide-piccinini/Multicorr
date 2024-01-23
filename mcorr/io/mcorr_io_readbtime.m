function [d,swapflag] = mcorr_io_readbtime(fid,forcebe)
% readbtime reads BTIME structure from current opened file and returns
%   D = [YEAR,DAY,HOUR,MINUTE,SECONDS]

Year        = fread(fid,1,'*uint16');
DayOfYear   = fread(fid,1,'*uint16');
Hours       = fread(fid,1,'uint8');
Minutes     = fread(fid,1,'uint8');
Seconds     = fread(fid,1,'uint8');
fseek(fid,1,0); % skip 1 byte (unused)
Seconds0001 = fread(fid,1,'*uint16');

% Automatic detection of little/big-endian encoding
% -- by F. Beauducel, March 2014 --
%
% If the 2-byte day is >= 512, the file is not opened in the correct
% endianness. If the day is 1 or 256, there is a possible byte-swap and we
% need to check also the year; but we need to consider what is a valid year:
% - years from 1801 to 2047 are OK (swapbytes >= 2312)
% - years from 2048 to 2055 are OK (swapbytes <= 1800)
% - year 2056 is ambiguous (swapbytes = 2056)
% - years from 2057 to 2311 are OK (swapbytes >= 2312)
% - year 1799 is ambiguous (swapbytes = 1799)
% - year 1800 is suspicious (swapbytes = 2055)
%
% Thus, the only cases for which we are 'sure' there is a byte-swap, are:
% - day >= 512
% - (day == 1 or day == 256) and (year < 1799 or year > 2311)
%
% Note: in IRIS libmseed, the test is only year>2050 or year<1920.
if ~forcebe && (DayOfYear >= 512 || (ismember(DayOfYear,[1,256]) && (Year > 2311 || Year < 1799)))
    swapflag = 1;
    Year = swapbytes(Year);
    DayOfYear = swapbytes(DayOfYear);
    Seconds0001 = swapbytes(Seconds0001);
else
    swapflag = 0;
end
d = [double(Year),double(DayOfYear),Hours,Minutes,Seconds + double(Seconds0001)/1e4];
