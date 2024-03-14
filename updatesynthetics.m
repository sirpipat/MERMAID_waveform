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
% Last modified by sirawich-at-princeton.edu, 03/14/2024

defval('model', 'ak135')

% read a synthetic SAC file
[SeisData, HdrData] = readsac(fname);

new_version = true;
% NEW VERSION
% Use TauP 1.3.0 or newer where we can specify receiver depth
if new_version
    tt = tauptime('mod', model, ...
                  'depth', HdrData.EVDP, ...
                  'gcarc', HdrData.GCARC, ...
                  'ph', 'p,s,P,S,PP,SS,PKP,SKS,PKIKP,SKIKS', ...
                  'stadepth', -HdrData.STEL/1000);
else
% OLD VERSION
% compute theoretical travel times at the ocean bottom below MERMAID.
% [lat lon] of the receiver is slightly shifted if incident angle is not
% close to zero.
    tt = taupPierce(model, HdrData.EVDP, ...
        'p,s,P,S,PP,SS,PKP,SKS,PKIKP,SKIKS', ...
        'sta', [HdrData.STLA HdrData.STLO], ...
        'evt', [HdrData.EVLA HdrData.EVLO], ...
        'pierce', -HdrData.STEL/1000, 'nodiscon');

    % remove all zero piercings
    for ii = 1:length(tt)
        index = length(tt(ii).pierce.p);
        while tt(ii).pierce.time(index) <= 0 && index > 1
            index = index - 1;
        end
        tt(ii).time = tt(ii).pierce.time(index);
        tt(ii).distance = tt(ii).pierce.distance(index);
    end
end

% keep only one arrival for each phase
ph = cell(size(tt));
for ii = 1:length(ph)
    if new_version
        tt(ii).phaseName = tt(ii).phase;
    end
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

if new_version
    HdrData.USER9 = tt(1).rayparameter;
else
    HdrData.USER9 = tt(1).rayParam;
end

% compute the time adjustment to account for receiver in Instaseis 
% is at the surface instead of the ocean bottom.
if strcmpi(model, 'ak135')
    vp = ak135('depths', -HdrData.STEL / 1000, 'dcbelow', false).vp;
elseif strcmpi(model, 'iasp91')
    vp = iasp91('depths', -HdrData.STEL / 1000, 'dcbelow', false).vp;
elseif strcmpi(model, 'prem')
    vp = prem('depths', -HdrData.STEL / 1000, 'dcbelow', false).vp;
else
    warning(['This model (%s) is not implemented for this function' ...
            'yet. PREM is used instead.\n'], model);
    vp = prem('depths', -HdrData.STEL / 1000, 'dcbelow', false).vp;
end
R_Earth = 6371;
theta_i = real(asin(HdrData.USER9 * vp / (R_Earth + HdrData.STEL / 1000)));
t_adjust = (HdrData.STEL / 1000) / (vp * cos(theta_i));

% adjust the reference time, begin time, and end time accordingly
% we use dt_B output here to acount for fractions of MSEC as well
[~, dt_ref_true] = gethdrinfo(HdrData);
dt_ref_true = dt_ref_true + seconds(t_adjust);
HdrData.NZYEAR = dt_ref_true.Year;
HdrData.NZJDAY = floor(days(dt_ref_true - datetime(dt_ref_true.Year, 1, ...
    0, 0, 0, 0, 0, 'TimeZone', 'UTC')));
HdrData.NZHOUR = dt_ref_true.Hour;
HdrData.NZMIN = dt_ref_true.Minute;
HdrData.NZSEC = floor(dt_ref_true.Second);
HdrData.NZMSEC = floor((dt_ref_true.Second - HdrData.NZSEC) * 1000);
B = dt_ref_true.Second - HdrData.NZSEC - (HdrData.NZMSEC / 1000);
E = HdrData.E + (B - HdrData.B);
HdrData.B = B;
HdrData.E = E;

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

% save SAC file
writesac(SeisData, HdrData, fname);
end