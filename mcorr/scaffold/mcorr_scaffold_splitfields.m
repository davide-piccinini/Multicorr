function c = mcorr_scaffold_splitfield(s,d)
% splitfield(S) splits string S of D-character separated field names

C = textscan(s,'%s','Delimiter',d);
c = C{1};
