function varargout = mcorr_io_rdmseed(varargin)
%RDMSEED Read miniSEED format file.
%   X = RDMSEED(F) reads file F and returns a M-by-1 structure X containing
%   M blocks ("data records") of a miniSEED file with headers, blockettes,
%   and data in dedicated fields, in particular, for each data block X(i):
%                t: time vector (DATENUM format)
%                d: data vector (double)
%       BLOCKETTES: existing blockettes (substructures)
%
%   Known blockettes are 100, 500, 1000, 1001 and 2000. Others will be
%   ignored with a warning message.
%
%   X = RDMSEED(F,ENCODINGFORMAT,WORDORDER,RECORDLENGTH), when file F does
%   not include the Blockette 1000 (like Seismic Handler outputs), specifies:
%       - ENCODINGFORMAT: FDSN code (see below); default is 10 = Steim-1;
%       - WORDORDER: 1 = big-endian (default), 0 = little-endian;
%       - RECORDLENGTH: must be a power of 2, at least 256 (default is 4096).
%   If the file contains Blockette 1000 (which is mandatory in the SEED
%   convention...), these 3 arguments are ignored except with 'force' option.
%
%   X = RDMSEED without input argument opens user interface to select the
%   file from disk.
%
%   [X,I] = RDMSEED(...) returns a N-by-1 structure I with N the detected
%   number of different channels, and the following fields:
%       ChannelFullName: channel name,
%           XBlockIndex: channel's vector index into X,
%            ClockDrift: vector of time interval errors, in seconds,
%                        between each data block (relative to sampling
%                        period). This can be compared to "Max Clock Drift"
%                        value of a Blockette 52.
%                           = 0 in perfect case
%                           < 0 tends to overlapping
%                           > 0 tends to gapping
%     OverlapBlockIndex: index of blocks (into X) having a significant
%                        overlap with previous block (less than 0.5
%                        sampling period).
%           OverlapTime: time vector of overlapped blocks (DATENUM format).
%         GapBlockIndex: index of blocks (into X) having a significant gap
%                        with next block (more than 0.5 sampling period).
%               GapTime: time vector of gapped blocks (DATENUM format).
%
%   RDMSEED(...) without output arguments plots the imported signal by
%   concatenating all the data records, in one single plot if single channel
%   is detected, or subplots for multiplexed file (limited to 10 channels).
%   Gaps are shown with red stars, overlaps with green circles.
%
%   [...] = RDMSEED(F,...,'be') forces big-endian reading (overwrites the
%   automatic detection of endianness coding, which fails in some cases).
%
%   [...] = RDMSEED(F,...,'notc') disables time correction.
%
%   [...] = RDMSEED(F,...,'nullhead') ignores null header (some files may
%   start with a series of null bytes).
%
%   [...] = RDMSEED(F,...,'plot') forces the plot with output arguments.
%
%   [...] = RDMSEED(F,...,'v') uses verbose mode (displays additional
%   information and warnings when necessary). Use 'vv' for extras, 'vvv'
%   for debuging.
%
%   Some instructions for usage of the returned structure:
%
%   - to get concatenated time and data vectors from a single-channel file:
%       X = rdmseed(f,'plot');
%       t = cat(1,X.t);
%       d = cat(1,X.d);
%
%   - to get the list of channels in a multiplexed file:
%       [X,I] = rdmseed(f);
%       char(I.ChannelFullName)
%
%   - to extract the station component n from a multiplexed file:
%       [X,I] = rdmseed(f);
%       k = I(n).XBlockIndex;
%       plot(cat(1,X(k).t),cat(1,X(k).d))
%       datetick('x')
%       title(I(n).ChannelFullName)
%
%   Known encoding formats are the following FDSN codes:
%        0: ASCII
%        1: 16-bit integer
%        2: 24-bit integer
%        3: 32-bit integer
%        4: IEEE float32
%        5: IEEE float64
%       10: Steim-1
%       11: Steim-2
%       12: GEOSCOPE 24-bit (untested)
%       13: GEOSCOPE 16/3-bit gain ranged
%       14: GEOSCOPE 16/4-bit gain ranged
%       19: Steim-3 (alpha and untested)
%
%   See also MKMSEED to export data in miniSEED format.
%
%
%   Author: François Beauducel <beauducel@ipgp.fr>
%       Institut de Physique du Globe de Paris
%   Created: 2010-09-17
%   Updated: 2018-08-09
%
%   Acknowledgments:
%       Ljupco Jordanovski, Jean-Marie Saurel, Mohamed Boubacar, Jonathan Berger,
%       Shahid Ullah, Wayne Crawford, Constanza Pardo, Sylvie Barbier,
%       Robert Chase, Arnaud Lemarchand, Alexandre Nercessian.
%
%       Special thanks to Martin Mityska who also inspired me with his ingenious
%       ReadMSEEDFast.m function.
%
%   References:
%       IRIS (2010), SEED Reference Manual: SEED Format Version 2.4, May 2010,
%         IFDSN/IRIS/USGS, http://www.iris.edu
%       Trabant C. (2010), libmseed: the Mini-SEED library, IRIS DMC.
%       Steim J.M. (1994), 'Steim' Compression, Quanterra Inc.

%   History:
%
%       [2018-08-09]
%           - MAJOR CODE UPDATE: now processes the binary data in memory
%             after a global file reading.
%           - removes all global variables.
%       [2017-11-21]
%           - adds option 'nullhead' to bypass null bytes header.
%       [2015-01-05]
%           - fixes a bug when a data block has 0 sample declared in header
%             but some data in the record (STEIM-1/2 coding).
%       [2014-06-29]
%           - 24-bit uncompressed format tested (bug correction), thanks to
%             Arnaud Lemarchand.
%       [2014-05-31]
%           - applies the time correction to StartTime and X.t (if needed).
%           - new option 'notc' to disable time correction.
%           - Geoscope 16/4 format passed real data archive tests.
%           - fixes a problem when plotting multiplexed channels (thanks to
%             Robert Chase).
%       [2014-03-14]
%           - Improved endianness automatic detection (see comments).
%           - Accepts mixed little/big endian encoding in a single file.
%           - minor fixes.
%       [2013-10-25]
%           - Due to obsolete syntax of bitcmp(0,N) in R2013b, replaces all
%             by: 2^N-1 (which is much faster...)
%       [2013-02-15]
%           - Tests also DayOfYear in header to determine automatically
%             little-endian coding of the file.
%           - Adds option 'be' to force big-endian reading (overwrites
%             automatic detection).
%       [2012-12-21]
%           - Adds a verbose mode
%       [2012-04-21]
%           - Correct bug with Steim + little-endian coding
%             (thanks to Shahid Ullah)
%       [2012-03-21]
%           - Adds IDs for warning messages
%       [2011-11-10]
%           - Correct bug with multiple channel name length (thanks to
%             Jonathan Berger)
%       [2011-10-27]
%           - Add LocationIdentifier to X.ChannelFullName
%       [2011-10-24]
%           - Validation of IEEE double encoding (with PQL)
%           - Import/plot data even with file integrity problem (like PQL)
%       [2011-07-21]
%           - Validation of ASCII encoding format (logs)
%           - Blockettes are now stored in substructures below a single
%             field X.BLOCKETTES
%           - Add import of blockettes 500 and 2000
%           - Accept multi-channel files with various data coding
%       [2010-10-16]
%           - Alpha-version of Steim-3 decoding...
%           - Extend output parameters with channel detection
%           - Add gaps and overlaps on plots
%           - Add possibility to force the plot
%       [2010-10-02]
%           - Add the input formats for GEOSCOPE multiplexed old data files
%           - Additional output argument with gap and overlap analysis
%           - Create a plot when no output argument are specified
%           - Optimize script coding (30 times faster STEIM decoding!)
%       [2010-09-28]
%           - Correction of a problem with STEIM-1 nibble 3 decoding (one
%             32-bit difference)
%           - Add reading of files without blockette 1000 with additional
%             input arguments (like Seismic Handler output files).
%           - Uses warning() function instead of fprintf().
%
%   Copyright (c) 2018, François Beauducel, covered by BSD License.
%   All rights reserved.
%
%   Redistribution and use in source and binary forms, with or without
%   modification, are permitted provided that the following conditions are
%   met:
%
%      * Redistributions of source code must retain the above copyright
%        notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in
%        the documentation and/or other materials provided with the distribution
%
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%   POSSIBILITY OF SUCH DAMAGE.

if nargin > 6
    error('Too many input arguments.')
end

% default input arguments
makeplot = 0;   % make plot flag
verbose = 0;    % verbose flag/level
forcebe = 0;    % force big-endian
ef = 10;        % encoding format default
wo = 1;         % word order default
rl = 2^12;      % record length default
force = 0;      % force input argument over blockette 1000 (UNDOCUMENTED)
notc = 0;       % force no time correction (over ActivityFlags)
nullhead = 0;   % allow null bytes before header

if nargin < 1
    [filename,pathname] = uigetfile('*','Please select a miniSEED file...');
    f = fullfile(pathname,filename);
else
    f = varargin{1};
end

if ~ischar(f) || ~exist(f,'file')
    error('File %s does not exist.',f);
end

if nargin > 1
    verbose = any(strcmpi(varargin,'v')) + 2*any(strcmpi(varargin,'vv')) ...
              + 3*any(strcmpi(varargin,'vvv'));
    makeplot = any(strcmpi(varargin,'plot'));
    forcebe = any(strcmpi(varargin,'be'));
    notc = any(strcmpi(varargin,'notc'));
    force = any(strcmpi(varargin,'force'));
    nullhead = any(strcmpi(varargin,'nullhead'));
end
nargs = (makeplot>0) + (verbose>0) + (forcebe>0) + (notc>0) + (force>0) ...
     + (nullhead>0);


if nargin > (1 + nargs)
    ef = varargin{2};
    if ~isnumeric(ef) || ~any(ef==[0:5,10:19,30:33])
        error('Argument ENCODINGFORMAT must be a valid FDSN code value.');
    end
end

if nargin > (2 + nargs)
    wo = varargin{3};
    if ~isnumeric(wo) || (wo ~= 0 && wo ~= 1)
        error('Argument WORDORDER must be 0 or 1.');
    end
end

if nargin > (3 + nargs)
    rl = varargin{4};
    if ~isnumeric(rl) || rl < 256 || rem(log(rl)/log(2),1) ~= 0
        error('Argument RECORDLENGTH must be a power of 2 and greater or equal to 256.');
    end
end

if nargout == 0
    makeplot = 1;
end

% sensible limits for multiplexed files
max_channels = 20;  % absolute max number of channels to plot
max_channel_label = 6;  % max. number of channels for y-labels

% file is opened in Big-Endian encoding (this is encouraged by SEED)
fid = fopen(f,'rb','ieee-be');
le = 0;
offset = 0;

% --- tests if the header is mini-SEED
% the 7th character must be one of the "data header/quality indicator", usually 'D'
header = fread(fid,20,'*char');
if ~ismember(header(7),'DRMQ')
    if ismember(header(7),'VAST')
        error('File seems to be a SEED Volume. Cannot read it.');
    else
        if header(1)==0
            if nullhead
                if verbose
                    fprintf('Null header option: bypassing...');
                end
                c = 0;
                fseek(fid,0,'bof');
                while c==0
                    c = fread(fid,1,'*char');
                    offset = offset + 1;
                end
                if verbose
                    fprintf(' %d null bytes.\n',offset);
                end
                header = fread(fid,6,'*char');
                if ~ismember(header(6),'DRMQ')
                    error('File is not in mini-SEED format. Cannot read it.');
                else
                    offset = offset - 1;
                end
            else
                error('File starts with null bytes... if you believe it is still a miniseed file, try the ''nullhead'' option.');
            end
        else
            error('File is not in mini-SEED format. Cannot read it.');
        end
    end
end

i = 1;

% --- main loop that reads data records until the end of the file
while offset >= 0
    [X(i),offset] = read_data_record(f,fid,offset,le,ef,wo,rl,forcebe,verbose,notc,force);
    i = i + 1;
end

fclose(fid);

if nargout > 0
    varargout{1} = X;
end

% --- analyses data
if makeplot || nargout > 1

    % test if the file is multiplexed or a single channel
    un = unique(cellstr(char(X.ChannelFullName)));
    nc = numel(un);
    for i = 1:nc
        k = find(strcmp(cellstr(char(X.ChannelFullName)),un{i}));
        I(i).ChannelFullName = X(k(1)).ChannelFullName;
        I(i).XBlockIndex = k;
        I(i).ClockDrift = ([diff(cat(1,X(k).RecordStartTimeMATLAB));NaN]*86400 - cat(1,X(k).NumberSamples)./cat(1,X(k).SampleRate))./cat(1,X(k).NumberSamples);
        I(i).OverlapBlockIndex = k(find(I(i).ClockDrift.*cat(1,X(k).NumberSamples).*cat(1,X(k).SampleRate) < -.5) + 1);
        I(i).OverlapTime = cat(1,X(I(i).OverlapBlockIndex).RecordStartTimeMATLAB);
        I(i).GapBlockIndex = k(find(I(i).ClockDrift.*cat(1,X(k).NumberSamples).*cat(1,X(k).SampleRate) > .5) + 1);
        I(i).GapTime = cat(1,X(I(i).GapBlockIndex).RecordStartTimeMATLAB);
    end
end
if nargout > 1
    varargout{2} = I;
end

% --- plots the data
if makeplot

    figure

    xlim = [min(cat(1,X.t)),max(cat(1,X.t))];

    % test if all data records have the same length
    rl = unique(cat(1,X.DataRecordSize));
    if numel(rl) == 1
        rl_text = sprintf('%d bytes',rl);
    else
        rl_text = sprintf('%d-%d bytes',min(rl),max(rl));
    end

    % test if all data records have the same sampling rate
    sr = unique(cat(1,X.SampleRate));
    if numel(sr) == 1
        sr_text = sprintf('%g Hz',sr);
    else
        sr_text = sprintf('%d # samp. rates',numel(sr));
    end

    % test if all data records have the same encoding format
    ef = unique(cellstr(cat(1,X.EncodingFormatName)));
    if numel(ef) == 1
        ef_text = sprintf('%s',ef{:});
    else
        ef_text = sprintf('%d different encod. formats',numel(ef));
    end

    if nc == 1
        plot(cat(1,X.t),cat(1,X.d))
        hold on
        for i = 1:length(I.GapBlockIndex)
            plot(I.GapTime(i),X(I.GapBlockIndex(i)).d(1),'*r')
        end
        for i = 1:length(I.OverlapBlockIndex)
            plot(I.OverlapTime(i),X(I.OverlapBlockIndex(i)).d(1),'og')
        end
        hold off
        set(gca,'XLim',xlim)
        datetick('x','keeplimits')
        grid on
        xlabel(sprintf('Time\n(%s to %s)',datestr(xlim(1)),datestr(xlim(2))))
        ylabel('Counts')
        title(sprintf('mini-SEED file "%s"\n%s (%d rec. @ %s - %g samp. @ %s - %s)', ...
            f,un{1},length(X),rl_text,numel(cat(1,X.d)),sr_text,ef_text),'Interpreter','none')
    else
        % plot is done only for real data channels...
        if nc > max_channels
            warning('Plot has been limited to %d channels (over %d). See help to manage multiplexed file.', ...
                max_channels,nc);
            nc = max_channels;
        end
        for i = 1:nc
            subplot(nc*2,1,i*2 + (-1:0))
            k = I(i).XBlockIndex;
            if ~any(strcmp('ASCII',cellstr(cat(1,X(k).EncodingFormatName))))
                plot(cat(1,X(k).t),cat(1,X(k).d))
                hold on
                for ii = 1:length(I(i).GapBlockIndex)
                    if ~isempty(X(I(i).GapBlockIndex(ii)).d)
                        plot(I(i).GapTime(ii),X(I(i).GapBlockIndex(ii)).d,'r')
                    else
                        plot(repmat(I(i).GapTime(ii),1,2),ylim,'r')
                    end
                end
                for ii = 1:length(I(i).OverlapBlockIndex)
                    if ~isempty(X(I(i).OverlapBlockIndex(ii)).d)
                        plot(I(i).OverlapTime(ii),X(I(i).OverlapBlockIndex(ii)).d,'g')
                    else
                        plot(repmat(I(i).OverlapTime(ii),1,2),ylim,'g')
                    end
                end
                hold off
            end
            set(gca,'XLim',xlim,'FontSize',8)
            h = ylabel(un{i},'Interpreter','none');
            if nc > max_channel_label
                set(gca,'YTick',[])
                set(h,'Rotation',0,'HorizontalAlignment','right','FontSize',8)
            end
            datetick('x','keeplimits')
            set(gca,'XTickLabel',[])
            grid on
            if i == 1
                title(sprintf('mini-SEED file "%s"\n%d channels (%d rec. @ %s - %g data - %s - %s)', ...
                    f,length(un),length(X),rl_text,numel(cat(1,X(k).d)),sr_text,ef_text),'Interpreter','none')
            end
            if i == nc
                datetick('x','keeplimits')
                xlabel(sprintf('Time\n(%s to %s)',datestr(xlim(1)),datestr(xlim(2))))
            end
        end
        v = version;
        if str2double(v(1))>=7
            linkaxes(findobj(gcf,'type','axes'),'x')
        end
    end
end
