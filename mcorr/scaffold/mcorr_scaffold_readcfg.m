function CFG=mcorr_scaffold_readcfg(inputfile)

% usage: CFG=mcorr_scaffold_readcfg(inputfile)
% Read configuration file and create a CFG structure
%
% Config file should be in the form 
% varname1 value1
% varname2 value2
%
% CFG structure will be:
% CFG.varname1.value1
% CFG.varname2.value2

% Davide Piccinini 5 Oct 2022

fprintf('Reading configuration from: %s\n',inputfile)
if exist (inputfile)~=2
    fprintf('No config file found\n');
    fprintf('** STOP **\n');
    CFG=struct;
    return
end

CFG=struct;
T = readtable(inputfile,'Delimiter','=', ...
              'EmptyColumnRule', 'skip', 'FileType', 'text');

[r,c] = size(T);
for ii = 1:r
    CFG=setfield(CFG,T.Var1{ii},T.Var2(ii));
end
