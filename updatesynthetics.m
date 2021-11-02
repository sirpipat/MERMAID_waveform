function updatesynthetics(fname, model)
% UPDATESYNTHETICS(fname, model)
%
% Takes a synthetic SAC file as an input. Then, computes theoretical
% travel times and ray parameters using TauP given the model name and
% source and station locations. Finally, saved SAC files with
% updated headers.
%
% INPUT:
% fname         full filename to a synthetic SAC file
% model         name of the model [default: 'ak135']
%
% It updates the following header fields:
% Tn            arrival time of n-th phase
% KTn           phase name of n-th phase
% USER9         ray parameter of the first arrival phase
%
% Last modified by sirawich-at-princeton.edu, 11/01/2021

defval('model', 'ak135')

% read a synthetic SAC file
[SeisData, HdrData] = readsac(fname);

% compute theoretical travel times
tt = taupTime(model, HdrData.EVDP, 'p,s,P,S,PP,SS,PcP,ScS,Pdiff,Sdiff,PKP,SKS,PKIKP,SKIKS', ...
    'sta', [HdrData.STLA HdrData.STLO], ...
    'evt', [HdrData.EVLA HdrData.EVLO]);

% keep only one arrival for each phase
ph = cell(size(tt));
for ii = 1:length(ph)
    ph{ii} = tt(ii).phaseName;
end
[~, ia] = unique(ph);
tt = tt(ia);

% sort the arrivals by time
tp = zeros(size(tt));
for ii = 1:length(tp)
    tp(ii) = tt(ii).time;
end
[~, is] = sort(tp);
tt = tt(is);

% clear the existing arrival-time tags
HdrData.T0 = -12345;
HdrData.T1 = -12345;
HdrData.T2 = -12345;
HdrData.T3 = -12345;
HdrData.T4 = -12345;
HdrData.T5 = -12345;
HdrData.T6 = -12345;
HdrData.T7 = -12345;
HdrData.T8 = -12345;
HdrData.T9 = -12345;

HdrData.KT0 = '-12345  ';
HdrData.KT1 = '-12345  ';
HdrData.KT2 = '-12345  ';
HdrData.KT3 = '-12345  ';
HdrData.KT4 = '-12345  ';
HdrData.KT5 = '-12345  ';
HdrData.KT6 = '-12345  ';
HdrData.KT7 = '-12345  ';
HdrData.KT8 = '-12345  ';
HdrData.KT9 = '-12345  ';

% update header
for ii = 1:length(tt)
    switch ii
        case 1
            HdrData.T0 = tt(ii).time;
            HdrData.KT0 = tt(ii).phaseName;
        case 2
            HdrData.T1 = tt(ii).time;
            HdrData.KT1 = tt(ii).phaseName;
        case 3
            HdrData.T2 = tt(ii).time;
            HdrData.KT2 = tt(ii).phaseName;
        case 4
            HdrData.T3 = tt(ii).time;
            HdrData.KT3 = tt(ii).phaseName;
        case 5
            HdrData.T4 = tt(ii).time;
            HdrData.KT4 = tt(ii).phaseName;
        case 6
            HdrData.T5 = tt(ii).time;
            HdrData.KT5 = tt(ii).phaseName;
        case 7
            HdrData.T6 = tt(ii).time;
            HdrData.KT6 = tt(ii).phaseName;
        case 8
            HdrData.T7 = tt(ii).time;
            HdrData.KT7 = tt(ii).phaseName;
        case 9
            HdrData.T8 = tt(ii).time;
            HdrData.KT8 = tt(ii).phaseName;
        case 10
            HdrData.T9 = tt(ii).time;
            HdrData.KT9 = tt(ii).phaseName;
        otherwise
            break
    end
end
HdrData.USER9 = tt(1).rayParam;

% save SAC file
writesac(SeisData, HdrData, fname);
end