% read in the data
clear all;clc;close all;
mon = [31,29,31,30,31,30,31,31,30,31,30,31];
mon_leap_sum = cumsum(mon);
mon = [31,28,31,30,31,30,31,31,30,31,30,31];
mon_sum = cumsum(mon);
delimiterIn = ',';
degree = sprintf('%c', char(176));
make_plot = 1;
% count the number of site years
count = 1;
Tair_ts = [];
Tair_d_ts = [];
Tsoil1_ts = [];
Tsoil1_d_ts = [];
Tsoil2_ts = [];
Tsoil2_d_ts = [];
FCH4_ts = [];
FCH4_d_ts = [];
nt_ts = [];
nt_d_ts = [];
hys_ts = [];
hys_no_IAV_ts = [];
aggregate = [];
hourly_data = [];

% variable names for analysis
var_wanted = ["TIMESTAMP_START","TA_F","TA_F_QC","P_F","P_F_QC","VPD_F",...
    "VPD_F_QC","NETRAD","NETRAD_QC","SW_IN_F_MDS","SW_IN_F_MDS_QC",...
    "WS_F","WS_F_QC","USTAR","USTAR_QC","PPFD_IN","PPFD_IN_QC",...
    "CO2_F_MDS","CO2_F_MDS_QC","TS_F_MDS_1","TS_F_MDS_1_QC",...
    "SWC_F_MDS_1","SWC_F_MDS_1_QC","G_F_MDS","G_F_MDS_QC","LE_F_MDS",...
    "LE_F_MDS_QC","H_F_MDS","H_F_MDS_QC","PA_F","PA_F_QC",...
    "NEE_VUT_REF","NEE_VUT_REF_QC",...
    "GPP_NT_VUT_REF","RECO_NT_VUT_REF"];
% read in the folders
path = '../site_level_measurements/fluxnet2015/';
dir=folderFiles([path],'FLX_*FLUXNET2015*');
% read in daily data
n_folder = length(dir(:,1));
site_included = [];
for site = 1:n_folder
    % locate the half hourly file
    file = folderFiles([path dir(site,:) '/'],...
        'FLX_*FLUXNET2015*FULLSET_HH*.csv');
    if (length(file)==0) % try hourly file
        file = folderFiles([path dir(site,:) '/'],...
            'FLX_*FLUXNET2015*FULLSET_HR*.csv');
    end
    % read in texts and values as a table
    T = readtable([path dir(site,:) '/' file]);
    % read in the data header
    fileID = fopen([path dir(site,:) '/' file]);
    tline = fgets(fileID);
    fclose(fileID);
    tline = strsplit(tline,',');
    header = string(erase(tline,'"'));
    hourly_tmp = [];
    for variables = 1:length(var_wanted) % loop through variables
        if (sum(contains(header,var_wanted(variables)))>=1) 
            % if the requested variable exists
            idx = find(strcmpi(header,var_wanted(variables)));
            tmp = T{:,idx};
            tmp(tmp==-9999) = nan;
            % variable will be stored as a cell array if it mixes texts &
            % numerics
            if (iscell(tmp))
                % find nan values in the array
                idx_NA=find(contains(tmp,'NA'));
                idx_full = logical(1:length(tmp))';
                idx_full(idx_NA) = 0;
                xxx = find(idx_full==1);
                yyy = find(idx_full==0);
                dim = size(tmp);
                data = zeros(dim);
                data(xxx) = cellfun(@str2num,tmp(xxx));
                data(yyy) = nan;
            else
                data = tmp;
            end
            hourly_tmp(1:length(T{:,idx}), variables) = data;
        else % the requested data is not available
            hourly_tmp(1:length(T{:,idx}), variables) = nan;
        end
            
    end
    % identify measurement length
    col = find(strcmpi("TIMESTAMP_START",var_wanted));
    tmp = floor(hourly_tmp(:,col)/(10^8));
    % only include data when there are >= 10 site-years
    if (length(unique(tmp))>=10)
        site_included = [site_included; site];
        hourly_tmp(:,length(var_wanted)+1) = site;
        hourly_data = [hourly_data; hourly_tmp];
    end
end

% read in site characteristics
path = '../site_level_measurements/fluxnet2015/';
% read in texts and values as a table
file = 'FLX_AA-Flx_BIF_LATEST.xlsx';
T = readtable([path file]);
clear IGBP LAT LONG
for site = 1:length(site_included)
    site_name = dir(site_included(site), 5:10);
    idx_site = find(strcmpi(T{:,1},site_name));
    tmp = T(idx_site,:);
    idx = find(strcmpi(tmp{:,4},'IGBP'));
    IGBP{site} = tmp{idx,5};
    idx = find(strcmpi(tmp{:,4},'LOCATION_LAT'));
    LAT(site) = mean(cellfun(@str2num,tmp{idx,5}));
    idx = find(strcmpi(tmp{:,4},'LOCATION_LONG'));
    LONG(site) = mean(cellfun(@str2num,tmp{idx,5}));
end
% convert IGBP to string array
IGBP = string(IGBP);
% convert LE to ET
var_test = ["TA_F"];
col_Tair = find(strcmpi(var_test(1),var_wanted)); % Tair
var_test = ["LE_F_MDS"];
col_LE = find(strcmpi(var_test(1),var_wanted)); % LE
% ET in kg H2O/m2/s
daily_ET = hourly_data(:,col_LE)/(.2491265*10^7);
% calculate saturation vapor pressure
daily_es = 6.11*10.^((7.5*hourly_data(:,col_Tair))./(237.3+hourly_data(:,col_Tair)));
% create time label
col = find(strcmpi("TIMESTAMP_START",var_wanted));
Time_Year = floor(hourly_data(:,col)/(10^8));
tmp = string(hourly_data(:,col));
clear Time_Month Time_Day Time_DOY Time_hr
for i = 1:length(tmp)
    mmdd = char(tmp(i));
    Time_Month(i) = str2num(mmdd(end-7:end-6));
    Time_Day(i) = str2num(mmdd(end-5:end-4));
    Time_hr(i) = str2num(mmdd(end-3:end));
    if(Time_Month(i)>1 & mod(Time_Year(i),4)~=0)
        Time_DOY(i) = mon_sum(Time_Month(i)-1)+Time_Day(i);
    elseif(Time_Month(i)>1 & mod(Time_Year(i),4)==0)
        Time_DOY(i) = mon_leap_sum(Time_Month(i)-1)+Time_Day(i);
    else
        Time_DOY(i) = Time_Day(i);
    end
end
Time = [];
Time = [Time_Year, Time_DOY', Time_hr'];

%% this section calculates WUE based on growing-season hourly measurements
% only include the measued data
clc;
close all;
% locate the variables: Tair; Precip; VPD; CO2 mole fraction; GPP; 
% Solar radiation; USTAR; NEE
var_test = ["TA_F","P_F","VPD_F",...
    "CO2_F_MDS","SW_IN_F_MDS","USTAR", "NEE_VUT_REF"];
var_idx = [];
for i = 1:length(var_test)
    var_idx = [var_idx; find(strcmpi(var_test(i),var_wanted))];
end
% measured data only
data = [];
for j = 1:length(var_idx)
    col = var_idx(j);
    data(:,j+1) = hourly_data(:,col);
    % QC flags, only use measured data
    idx = find(hourly_data(:,col+1)>0);
    data(idx,j+1) = nan;
end
% set time label
data(:,1) = Time(:,2);%Time_DOY;
var_test = ["RECO_NT_VUT_REF"];
var_idx = [];
for i = 1:length(var_test)
    var_idx = [var_idx; find(strcmpi(var_test(i),var_wanted))];
end
reco = hourly_data(:,var_idx(1));
% umolCO2 m-2 s-1 to gCm-2s-1
data(:,end+1) = reco*12*(10^-6);

% locate GPP
% umolCO2 m-2 s-1 for hourly data
var_test = ["GPP_NT_VUT_REF"];
var_idx = [];
for i = 1:length(var_test)
    var_idx = [var_idx; find(strcmpi(var_test(i),var_wanted))];
end
% GPP QC
gpp = hourly_data(:,var_idx(1));
% set GPP, umolCO2 m-2 s-1 to gCm-2s-1
data(:,end+1) = gpp*12*(10^-6);
% set ET
data(:,end+1) = daily_ET;
% include site-year info
data(:,end+1) = Time(:,1);%Time_Year;
% include site-hour info
data(:,end+1) = Time(:,3);%Time_hr;
% include site info
data(:,end+1) = hourly_data(:,end);
vpd_m = [];
Tair_m = [];
cco2_m = [];
% define daytime as solar radiation > 100 W/m2
sw_threshold = 100;
et_threshold = 10^-6;
gpp_threshold = 10^-6;
tair_threshold = 0;
ustar_threshold = 0.1;
wue_hourly_agg = [];
wue_daily_agg = [];
wue_monthly_agg = [];
grs_mean_agg_all = [];
gpp_pct = [];
et_pct = [];
vpd_pct = [];
% declare month DOY
for month = 2:12
    month_start(month) = mon_sum(month-1)+1;
    month_end(month) = mon_sum(month);
end
month_start(1) = 1;
month_end(1) = mon_sum(1);
pct_window = 5:5:95;
clear uwue iwue obs_year pct_value
% calculate WUE based on seasonal median/sums (1 WUE per growing season)
for site = 1:length(site_included)
    lat_tmp = LAT(site);
    if (IGBP(site)~='CRO')
        count = 1;
        idx = find(data(:,14)==site_included(site));
        data_site = data(idx,:);
        yr_range = [data_site(1, 12):data_site(end, 12)];
        if (length(yr_range)>0)
            yr_start(site) = yr_range(1);
            yr_end(site) = yr_range(end);
            for yr = 1:length(yr_range)
                idx = find(data_site(:,12)==yr_range(yr));
                data_site_year = data_site(idx,:);
                % precip events
                precip = sum(data_site_year(:,3),'omitnan'); % total precip [mm]
                % precip events
                idx = find(data_site_year(:,3)>0);
                doy_tmp = unique(data_site_year(idx,1));
                % filter out data during and one day after a precip event
                for dd = 1:length(doy_tmp)
                    idx = find(data_site_year(:,1)==doy_tmp(dd)|...
                        data_site_year(:,1)==(doy_tmp(dd)+1));
                    data_site_year(idx,:) = nan;
                end
                % filter out low ustar
                idx = find(data_site_year(:,7)<ustar_threshold);
                data_site_year(idx,:) = nan;
                % filter out nighttime data
                idx = find(data_site_year(:,6)<sw_threshold);
                data_site_year(idx,:) = nan;
                % WUE, approach#2
                % condensation events
                idx = find(data_site_year(:,11)<et_threshold);
                data_site_year(idx,:) = nan;
                % set growing season as the period when daily GPP>5%
                % daily GPP max and daily Tair>0
                idx = isfinite(data_site_year(:, 1));
                days = unique(data_site_year(idx, 1));
                tmp = [];
                for day = 1:length(days)
                    idx = find(data_site_year(:, 1)==days(day));
                    if (length(idx)>0)
                        % DOY, Tair, GPP
                        xxx = [data_site_year(idx(1), 1), ...
                            mean(data_site_year(idx,2), 'omitnan'),...
                            sum(data_site_year(idx,10), 'omitnan')];
                        tmp = [tmp; xxx];
                    end
                end
                idx = find(tmp(:,2)<=0 | tmp(:,3)<max(tmp(:,3))*0.05);
                if (length(idx)>1)
                    days = unique(tmp(idx, 1));
                    for day = 1:length(days)
                        idx = find(data_site_year(:, 1)==days(day));
                        if (length(idx)>0)
                            data_site_year(idx,:) = nan;
                        end
                    end
                end
                % When GPP<0
                idx = find(data_site_year(:,10)<gpp_threshold);
                data_site_year(idx,:) = nan;
                % When Tair<0
                idx = find(data_site_year(:,2)<tair_threshold);
                data_site_year(idx,:) = nan;
                iwue_hourly = data_site_year(:,10).*data_site_year(:,4)./data_site_year(:,11);
                uwue_hourly = data_site_year(:,10).*(data_site_year(:,4).^.5)./data_site_year(:,11);
                % daily aggregate
                % output year, DOY, WUE, site, GPP, ET, VPD, Tair, [CO2],
                % USTAR, SW_IN, NEE,RECO
                idx = isfinite(data_site_year(:, 1));
                days = unique(data_site_year(idx, 1));
                tmp = [];
                for day = 1:length(days)
                    idx = find(data_site_year(:, 1)==days(day));
                    if (length(idx)>0)
                        xxx = [data_site_year(idx(1),12), data_site_year(idx(1), 1), ...
                        mean(uwue_hourly(idx), 'omitnan'), data_site_year(idx(1),14), ...
                        sum(data_site_year(idx,10), 'omitnan'), sum(data_site_year(idx,11), 'omitnan'), ...
                        mean(data_site_year(idx,4), 'omitnan'), mean(data_site_year(idx,2), 'omitnan'), ...
                        mean(data_site_year(idx,5), 'omitnan'), mean(data_site_year(idx,7), 'omitnan'), ...
                        mean(data_site_year(idx,6), 'omitnan'), sum(data_site_year(idx,8), 'omitnan'), ...
                        sum(data_site_year(idx,9), 'omitnan')];
                        xxx(3) = xxx(5)*(xxx(7)^0.5)/xxx(6); % uWUE from daily agg values
                        xxx(end+1) = xxx(5)*(xxx(7))/xxx(6); % WUEei from daily agg values
                        tmp = [tmp; xxx];
                    end
                end
                wue_daily = tmp;
                if (length(wue_daily)>0) % data exist
                    tmp = [data_site_year(:,12), data_site_year(:, 1), ...
                        uwue_hourly, data_site_year(:,14), data_site_year(:,10), ...
                        data_site_year(:,11), data_site_year(:,4), ...
                        data_site_year(:,2), data_site_year(:,5), ...
                        data_site_year(:,7), data_site_year(:,13), ...
                        data_site_year(:,6), data_site_year(:,8), ...
                        data_site_year(:,9), iwue_hourly];
                    wue_hourly_agg = [wue_hourly_agg; tmp];
                    wue_daily_agg = [wue_daily_agg; wue_daily];
                    % monthly aggregate
                    % output year, month, WUE, site, GPP, ET, VPD,
                    % Tair, [CO2], USTAR, SW_IN, NEE, RECO
                    tmp = [];
                    for month = 1:12
                        idx = find((data_site_year(:, 1)>=month_start(month))&...
                            (data_site_year(:, 1)<=month_end(month)));
                        if (length(idx)>0)
                            xxx = [data_site_year(idx(1),12), month, ...
                            mean(uwue_hourly(idx), 'omitnan'), data_site_year(idx(1),14), ...
                            sum(data_site_year(idx,10), 'omitnan'), sum(data_site_year(idx,11), 'omitnan'), ...
                            mean(data_site_year(idx,4), 'omitnan'), mean(data_site_year(idx,2), 'omitnan'), ...
                            mean(data_site_year(idx,5), 'omitnan'), mean(data_site_year(idx,7), 'omitnan'), ...
                            mean(data_site_year(idx,6), 'omitnan'), sum(data_site_year(idx,8), 'omitnan'),...
                            sum(data_site_year(idx,9), 'omitnan')];
                            xxx(end+1) = xxx(5)*(xxx(7))/xxx(6); % WUEei from monthly agg values
                            tmp = [tmp; xxx];
                        end
                    end
                    wue_monthly_agg = [wue_monthly_agg; tmp];
                    % output growing season mean values
                    tmp = data_site_year(:, [1:12, 14]);
                    grs_mean = mean(tmp, 'omitnan');
                    grs_mean(3) = precip;
                    grs_mean(14) = sum(data_site_year(:,10), 'omitnan'); % gpp sum
                    grs_mean(15) = sum(data_site_year(:,11), 'omitnan'); % ET sum
                    grs_mean(16) = median(uwue_hourly, 'omitnan'); % median WUE from hourly WUE
                    grs_mean(17) = sum(data_site_year(:,10), 'omitnan')...
                        *(mean(data_site_year(:,4), 'omitnan')^.5)/sum(data_site_year(:,11), 'omitnan'); % GRS WUE
                    grs_mean(18) = median(iwue_hourly, 'omitnan'); % median iWUE from hourly iWUE
                    grs_mean(19) = sum(data_site_year(:,10), 'omitnan')...
                        *mean(data_site_year(:,4), 'omitnan')/sum(data_site_year(:,11), 'omitnan'); % GRS WUE
                    wue = uwue_hourly;
                    for window = 1:length(pct_window)
                        % for each percentile
                        grs_mean(17+window) = prctile(wue,pct_window(window));
                    end
                    grs_mean_agg_all = [grs_mean_agg_all; grs_mean];
                end
            end
        end
    end
end

%% this section calculates WUE based on hourly measurements during summer
% only include the measued data
clc;
close all;
% locate the variables: Tair; Precip; VPD; CO2 mole fraction; GPP; 
% Solar radiation; USTAR; NEE
var_test = ["TA_F","P_F","VPD_F",...
    "CO2_F_MDS","SW_IN_F_MDS","USTAR", "NEE_VUT_REF"];
var_idx = [];
for i = 1:length(var_test)
    var_idx = [var_idx; find(strcmpi(var_test(i),var_wanted))];
end
% measured data only
data = [];
for j = 1:length(var_idx)
    col = var_idx(j);
    data(:,j+1) = hourly_data(:,col);
    % QC flags, only use measured and good quality gapfill data
    idx = find(hourly_data(:,col+1)>0);
    data(idx,j+1) = nan;
end
% set time label
data(:,1) = Time(:,2);%Time_DOY;
% locate RECO
% umolCO2 m-2 s-1 for hourly data
var_test = ["RECO_NT_VUT_REF"];
var_idx = [];
for i = 1:length(var_test)
    var_idx = [var_idx; find(strcmpi(var_test(i),var_wanted))];
end
reco = hourly_data(:,var_idx(1));
% umolCO2 m-2 s-1 to gCm-2s-1
data(:,end+1) = reco*12*(10^-6);

% locate GPP
% umolCO2 m-2 s-1 for hourly data
var_test = ["GPP_NT_VUT_REF"];
var_idx = [];
for i = 1:length(var_test)
    var_idx = [var_idx; find(strcmpi(var_test(i),var_wanted))];
end
% GPP QC
gpp = hourly_data(:,var_idx(1));
% set GPP, umolCO2 m-2 s-1 to gCm-2s-1
data(:,end+1) = gpp*12*(10^-6);
% set ET
data(:,end+1) = daily_ET;
% include site-year info
data(:,end+1) = Time(:,1);%Time_Year;
% include site-hour info
data(:,end+1) = Time(:,3);%Time_hr;
% include site info
data(:,end+1) = hourly_data(:,end);
vpd_m = [];
Tair_m = [];
cco2_m = [];
% define daytime as solar radiation > 100 W/m2
sw_threshold = 100;
et_threshold = 10^-6;
gpp_threshold = 10^-6;
tair_threshold = 0;
ustar_threshold = 0.1;
wue_hourly_summer = [];
wue_daily_summer = [];
wue_monthly_summer = [];
% declare month DOY
for month = 2:12
    month_start(month) = mon_sum(month-1)+1;
    month_end(month) = mon_sum(month);
end
month_start(1) = 1;
month_end(1) = mon_sum(1);
pct_window = 5:5:95;
clear uwue iwue obs_year pct_value
for site = 1:length(site_included)
    lat_tmp = LAT(site);
    % exclude croplands
    if (IGBP(site)~='CRO')
        count = 1;
        idx = find(data(:,14)==site_included(site));
        data_site = data(idx,:);
        yr_range = [data_site(1, 12):data_site(end, 12)];
        if (length(yr_range)>0)
            yr_start(site) = yr_range(1);
            yr_end(site) = yr_range(end);
            for yr = 1:length(yr_range)
                idx = find(data_site(:,12)==yr_range(yr));
                data_site_year = data_site(idx,:);
                % precip events
                precip = sum(data_site_year(:,3),'omitnan'); % total precip [mm]
                % precip events
                idx = find(data_site_year(:,3)>0);
                doy_tmp = unique(data_site_year(idx,1));
                % filter out data during and one day after a precip event
                for dd = 1:length(doy_tmp)
                    idx = find(data_site_year(:,1)==doy_tmp(dd)|...
                        data_site_year(:,1)==(doy_tmp(dd)+1));
                    data_site_year(idx,:) = nan;
                end
                % filter out low ustar
                idx = find(data_site_year(:,7)<ustar_threshold);
                data_site_year(idx,:) = nan;
                % filter out nighttime data
                idx = find(data_site_year(:,6)<sw_threshold);
                data_site_year(idx,:) = nan;
                % WUE, approach#2
                % condensation events
                idx = find(data_site_year(:,11)<et_threshold);
                data_site_year(idx,:) = nan;
                % define summer season
                if (lat_tmp>0)
                    summer_start = mon_sum(5)+1;
                    summer_end = mon_sum(8);
                    idx_other = find(data_site_year(:, 1)<summer_start | ...
                        data_site_year(:, 1)>summer_end);
                else
                    summer_start = mon_sum(11)+1;
                    summer_end = mon_sum(2);
                    idx_other = find(data_site_year(:, 1)<summer_start & ...
                        data_site_year(:, 1)>summer_end);
                end
                data_site_year(idx_other,:) = nan;
                % set growing season as the period when daily GPP>5%
                % daily GPP max and daily Tair>0
                idx = isfinite(data_site_year(:, 1));
                days = unique(data_site_year(idx, 1));
                tmp = [];
                for day = 1:length(days)
                    idx = find(data_site_year(:, 1)==days(day));
                    if (length(idx)>0)
                        % DOY, Tair, GPP
                        xxx = [data_site_year(idx(1), 1), ...
                            mean(data_site_year(idx,2), 'omitnan'),...
                            sum(data_site_year(idx,10), 'omitnan')];
                        tmp = [tmp; xxx];
                    end
                end
                idx = find(tmp(:,2)<=0 | tmp(:,3)<max(tmp(:,3))*0.05);
                if (length(idx)>1)
                    days = unique(tmp(idx, 1));
                    for day = 1:length(days)
                        idx = find(data_site_year(:, 1)==days(day));
                        if (length(idx)>0)
                            data_site_year(idx,:) = nan;
                        end
                    end
                end
                % When GPP<0
                idx = find(data_site_year(:,10)<gpp_threshold);
                data_site_year(idx,:) = nan;
                % When Tair<0
                idx = find(data_site_year(:,2)<tair_threshold);
                data_site_year(idx,:) = nan;
                iwue_hourly = data_site_year(:,10).*data_site_year(:,4)./data_site_year(:,11);
                uwue_hourly = data_site_year(:,10).*(data_site_year(:,4).^.5)./data_site_year(:,11);
                % daily aggregate
                % output year, DOY, WUE, site, GPP, ET, VPD, Tair, [CO2],
                % USTAR, SW_IN, NEE,RECO
                idx = isfinite(data_site_year(:, 1));
                days = unique(data_site_year(idx, 1));
                tmp = [];
                for day = 1:length(days)
                    idx = find(data_site_year(:, 1)==days(day));
                    if (length(idx)>0)
                        xxx = [data_site_year(idx(1),12), data_site_year(idx(1), 1), ...
                        mean(uwue_hourly(idx), 'omitnan'), data_site_year(idx(1),14), ...
                        sum(data_site_year(idx,10), 'omitnan'), sum(data_site_year(idx,11), 'omitnan'), ...
                        mean(data_site_year(idx,4), 'omitnan'), mean(data_site_year(idx,2), 'omitnan'), ...
                        mean(data_site_year(idx,5), 'omitnan'), mean(data_site_year(idx,7), 'omitnan'), ...
                        mean(data_site_year(idx,6), 'omitnan'), sum(data_site_year(idx,8), 'omitnan'), ...
                        sum(data_site_year(idx,9), 'omitnan')];
                        xxx(3) = xxx(5)*(xxx(7)^0.5)/xxx(6); % uWUE from daily agg values
                        xxx(end+1) = xxx(5)*(xxx(7))/xxx(6); % WUEei from daily agg values
                        tmp = [tmp; xxx];
                    end
                end
                wue_daily = tmp;
                if (length(wue_daily)>0) % data exist
                    tmp = [data_site_year(:,12), data_site_year(:, 1), ...
                        uwue_hourly, data_site_year(:,14), data_site_year(:,10), ...
                        data_site_year(:,11), data_site_year(:,4), ...
                        data_site_year(:,2), data_site_year(:,5), ...
                        data_site_year(:,7), data_site_year(:,13), ...
                        data_site_year(:,6), data_site_year(:,8), ...
                        data_site_year(:,9), iwue_hourly];
                    wue_hourly_summer = [wue_hourly_summer; tmp];
                    wue_daily_summer = [wue_daily_summer; wue_daily];
                    % monthly aggregate
                    % output year, month, WUE, site, GPP, ET, VPD,
                    % Tair, [CO2], USTAR, SW_IN, NEE, RECO
                    tmp = [];
                    for month = 1:12
                        idx = find((data_site_year(:, 1)>=month_start(month))&...
                            (data_site_year(:, 1)<=month_end(month)));
                        if (length(idx)>0)
                            xxx = [data_site_year(idx(1),12), month, ...
                            mean(uwue_hourly(idx), 'omitnan'), data_site_year(idx(1),14), ...
                            sum(data_site_year(idx,10), 'omitnan'), sum(data_site_year(idx,11), 'omitnan'), ...
                            mean(data_site_year(idx,4), 'omitnan'), mean(data_site_year(idx,2), 'omitnan'), ...
                            mean(data_site_year(idx,5), 'omitnan'), mean(data_site_year(idx,7), 'omitnan'), ...
                            mean(data_site_year(idx,6), 'omitnan'), sum(data_site_year(idx,8), 'omitnan'),...
                            sum(data_site_year(idx,9), 'omitnan')];
                            xxx(end+1) = xxx(5)*(xxx(7))/xxx(6); % WUEei from monthly agg values
                            tmp = [tmp; xxx];
                        end
                    end
                    wue_monthly_summer = [wue_monthly_summer; tmp];
                end
            end
        end
    end
end


%% analyze results for savana, forests and grasslands
clc;
IGBP_included = unique(IGBP);% % exclude croplands
idx = find(IGBP_included=='CRO'|IGBP_included=='WET');
IGBP_included(idx) = [];
grs_mean_subset = [];
LAT_subset = [];
LONG_subset = [];
IGBP_subset = [];
site_year = [];
for i = 1:length(IGBP_included)
    % find sites within an ecosystem type
    idx_IGBP = find(IGBP==IGBP_included(i));
    site_count(i) = length(idx_IGBP);
    for j = 1:length(idx_IGBP)
        site = site_included(idx_IGBP(j));
        % extract site data
        idx = find(grs_mean_agg_all(:,13)==site);        
        % include the site data if there's more than 10 site-years
        if (sum(isfinite(grs_mean_agg_all(idx,16)))>=10)
            IGBP_subset = [IGBP_subset; i];
            LAT_subset = [LAT_subset; LAT(idx_IGBP(j))];
            LONG_subset = [LONG_subset; LONG(idx_IGBP(j))];
            site_year = [site_year; length(idx)];
            tmp = [grs_mean_agg_all(idx,:), ones(length(idx),1)*i]; % attach type info
            grs_mean_subset = [grs_mean_subset; tmp];
        end
    end
end

%% plot the site map
clc;
close all
ax = worldmap('World');
load coastlines
geoshow(ax, coastlat,coastlon, 'DisplayType', 'polygon');
zoom(1.2)
idx_IBBP = unique(IGBP_subset);
h_type = [];
for i = 1:length(idx_IBBP)
    idx = find(IGBP_subset==idx_IBBP(i));
    for j = 1:length(idx)
        % site location
        ppp = plotm(LAT_subset(idx(j)), LONG_subset(idx(j)),'o','MarkerEdgeColor',...
            mcolor(i+2,:), 'MarkerFaceColor', mcolor(i+2,:), 'Markersize',5);
        if (j==1)
            h_type = [h_type; ppp];
        end
        hold on
        % site-year info
        ppp = plotm(LAT_subset(idx(j)), LONG_subset(idx(j)),'o','MarkerEdgeColor',...
            'k', 'MarkerFaceColor', 'none', 'Markersize',max(5, site_year(idx(j))));
        mlabel off;
        plabel off;
    end
end
tmp = IGBP_included(idx_IBBP);
lgd = legend(h_type,...
            tmp,...
            'orientation','horizontal', 'box','off');
lgd.FontSize = 12;
lgd.Position = [0.3 0.1 0.3643 0.0452];

figname=['../plots_hourly/uWUE_trend_full_v7/site_map.jpeg'];
print('-dpng','-r300',[figname]);

%% compare WUE_median and WUE_aggregated during growing season and summer
for season = 1:2
    if (season==1)
        wue_data = wue_hourly_agg;
        season_name = 'grs';
    elseif (season==2)
        wue_data = wue_hourly_summer;
        season_name = 'sum';
    end
    for wue_type = 1:2
        if (wue_type==1)
            yrange = [0 30];
            wue = 'uwue';
            col_wue = 3;
            wue_name = '{\it uWUE}';
            units = '(gC hPa^{0.5} kg H_2O^{-1})';
            wue_med_lgd = '{\it uWUE_{median}}';
            wue_agg_lgd = '{\it uWUE_{aggregate}}';
        elseif (wue_type==2)
            if (season==1)
                yrange = [0 80];
            elseif (season==2)
                yrange = [0 100];
            end
            wue = 'iwue';
            col_wue = 15;
            wue_name = '{\it WUE_{ei}}';
            units = '(gC hPa kg H_2O^{-1})';
            wue_med_lgd = '{\it WUE_{ei, median}}';
            wue_agg_lgd = '{\it WUE_{ei, aggregated}}';
        end
        clc;
        close all
        clear comp_wue_mean comp_wue_agg comp_sen_wue_mean comp_sen_wue_agg IGBP_name
        clear comp_std_wue_mean comp_std_wue_agg comp_sen_wue_up_mean comp_sen_wue_up_agg
        clear comp_sen_wue_low_mean comp_sen_wue_low_agg comp_data_agg comp_data_mean
        clear comp_sen_wue_mean_raw comp_sen_wue_agg_raw IGBP_name_w_yr
        clear comp_sen_wue_up_mean_raw comp_sen_wue_up_agg_raw
        clear comp_sen_wue_low_mean_raw comp_sen_wue_low_agg_raw
        clear comp_site_years comp_siteID
        id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", ...
            "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", ...
            "(q)", "(r)", "(s)", "(t)", "(u)", "(v)", "(w)", "(x)"];
        % analyze the top X% windows
        pct_window = [95];
        id_trend = [];
        IGBP_included = unique(IGBP);% % exclude croplands
        idx = find(IGBP_included=='CRO'|IGBP_included=='WET');
        IGBP_included(idx) = [];
        data_mean_agg = [];
        data_agg_agg = [];
        h_agg = [];
        sig_agg = [];
        for i = 1:length(IGBP_included)
            % find sites within an ecosystem type
            idx_IGBP = find(IGBP==IGBP_included(i));
            % for sites under each IGBP
            clear site_list
            annual_change = [];
            plot_site = [];
            wue_mean_agg = [];
            wue_agg_agg = [];
            sen_wue_mean_agg = [];
            sen_wue_agg_agg = [];
            sen_wue_mean_raw_agg = [];
            sen_wue_agg_raw_agg = [];
            std_wue_mean_agg = [];
            std_wue_agg_agg = [];
            sen_wue_up_mean = [];
            sen_wue_low_mean = [];
            sen_wue_up_agg = [];
            sen_wue_low_agg = [];
            sen_wue_up_mean_raw = [];
            sen_wue_low_mean_raw = [];
            sen_wue_up_agg_raw = [];
            sen_wue_low_agg_raw = [];
            pct_wue_mean_agg = [];
            pct_wue_agg_agg = [];
            yr_count = 0;
            for j = 1:length(idx_IGBP) % loop through sites
                clf;
                plot_year_idx = [];
                wue_mean = [];
                wue_median = [];
                wue_agg = [];
                sen_wue_mean = [];
                sen_wue_agg = [];
                sen_wue_mean_raw = [];
                sen_wue_agg_raw = [];
                valid_yr = [];
                siteID = dir(site_included(idx_IGBP(j)), 5:10);
                idx = find(wue_data(:,4)==site_included(idx_IGBP(j)));
                data_tmp = wue_data(idx, :);
                % define growing season as GPP>0 & Tair>0
                idx = find(data_tmp(:, 5)>0 & data_tmp(:, 8)>0);
                grs_data = data_tmp(idx,:);
                yr_range = unique(grs_data(:, 1));
                if (length(yr_range)>=10)
                    yr_count = yr_count+length(yr_range);
                    plot_site = [plot_site; j];
                    IGBP_name{i} = char(IGBP_included(i));
                    clear grs_sum
                    % for individual years
                    for yr = 1:length(yr_range) % loop through years
                        idx = find(grs_data(:, 1)==yr_range(yr));
                        data_site_year = grs_data(idx, :);
                        wue_mean = [wue_mean; mean(data_site_year(:, col_wue),'omitnan')];
                        wue_median = [wue_median; median(data_site_year(:, col_wue),'omitnan')];
                        % GRS WUE from cumulated GPP and ET
                        if (wue_type==1)
                            tmp = sum(data_site_year(:, 5),'omitnan')*(mean(data_site_year(:, 7),'omitnan')^.5)...
                                /sum(data_site_year(:, 6),'omitnan');
                        elseif (wue_type==2)
                            tmp = sum(data_site_year(:, 5),'omitnan')*(mean(data_site_year(:, 7),'omitnan'))...
                                /sum(data_site_year(:, 6),'omitnan');
                        end
                        wue_agg = [wue_agg; tmp];
                        if (sum(isfinite(data_site_year(:, col_wue)))>0)
                            valid_yr = [valid_yr, yr];
                        end
                    end
                    % temporal trend with Sen's method, normalized
                    % mean
                    if (length(valid_yr)>=10)
                        datain = nan(length(valid_yr),2);
                        datain(:,1) = yr_range(valid_yr);
                        datain(:,2) = (wue_median(valid_yr)-mean(wue_median(valid_yr)))/mean(wue_median(valid_yr));
                        [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] =...
                                ktaub(datain, 0.05);
                        h_agg = [h_agg; h];
                        sig_agg = [sig_agg; sig];
                        sen_wue_mean = [sen_wue_mean; sen];
                        sen_wue_up_mean = [sen_wue_up_mean; CIupper];
                        sen_wue_low_mean = [sen_wue_low_mean; CIlower];
                        % variability among years at the same site
                        std_wue_mean_agg = [std_wue_mean_agg; std(wue_mean(valid_yr))];
                        % pct change per year
                        tmp = diff(wue_median(valid_yr))./diff(yr_range(valid_yr))/mean(wue_median(valid_yr))*100;
                        pct_wue_mean_agg = [pct_wue_mean_agg; tmp];
                        datain = nan(length(valid_yr),2);
                        datain(:,1) = yr_range(valid_yr);
                        datain(:,2) = (wue_agg(valid_yr)-mean(wue_agg(valid_yr)))/mean(wue_agg(valid_yr));
                        [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] =...
                                ktaub(datain, 0.05);
                        sen_wue_agg = [sen_wue_agg; sen];
                        sen_wue_up_agg = [sen_wue_up_agg; CIupper];
                        sen_wue_low_agg = [sen_wue_low_agg; CIlower];
                        std_wue_agg_agg = [std_wue_agg_agg; std(wue_agg(valid_yr))];
                        % pct change per year
                        tmp = diff(wue_agg(valid_yr))./diff(yr_range(valid_yr))/mean(wue_agg(valid_yr))*100;
                        pct_wue_agg_agg = [pct_wue_agg_agg; tmp];
                        % temporal trend with Sen's method
                        datain = nan(length(valid_yr),2);
                        datain(:,1) = yr_range(valid_yr);
                        datain(:,2) = wue_mean(valid_yr);
                        data_mean_agg = [data_mean_agg; datain];
                        [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] =...
                                ktaub(datain, 0.05);
                        sen_wue_mean_raw = [sen_wue_mean_raw; sen];
                        sen_wue_up_mean_raw = [sen_wue_up_mean_raw; CIupper];
                        sen_wue_low_mean_raw = [sen_wue_low_mean_raw; CIlower];
                        % agg
                        datain = nan(length(valid_yr),2);
                        datain(:,1) = yr_range(valid_yr);
                        datain(:,2) = wue_agg(valid_yr);
                        data_agg_agg = [data_agg_agg; datain];
                        [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] =...
                                ktaub(datain, 0.05);
                        sen_wue_agg_raw = [sen_wue_agg_raw; sen];
                        sen_wue_up_agg_raw = [sen_wue_up_agg_raw; CIupper];
                        sen_wue_low_agg_raw = [sen_wue_low_agg_raw; CIlower];
                        wue_mean_agg = [wue_mean_agg; wue_median];
                        wue_agg_agg = [wue_agg_agg; wue_agg];
                        sen_wue_mean_agg = [sen_wue_mean_agg; sen_wue_mean];
                        sen_wue_agg_agg = [sen_wue_agg_agg; sen_wue_agg];
                        sen_wue_mean_raw_agg = [sen_wue_mean_raw_agg; sen_wue_mean_raw];
                        sen_wue_agg_raw_agg = [sen_wue_agg_raw_agg; sen_wue_agg_raw];
                        comp_site_years(j,i) = length(valid_yr);
                        comp_siteID{j,i} = siteID;
                    end
                end
            end
            IGBP_name_w_yr{i} = [char(IGBP_included(i)) ' (' int2str(sum(isfinite(wue_mean_agg))) ')'];
            % store in type array
            comp_wue_mean(1:length(wue_mean_agg),i) = wue_mean_agg;
            comp_wue_agg(1:length(wue_agg_agg),i) = wue_agg_agg;
            comp_sen_wue_mean(1:length(sen_wue_mean_agg),i) = sen_wue_mean_agg;
            comp_sen_wue_agg(1:length(sen_wue_agg_agg),i) = sen_wue_agg_agg;
            comp_sen_wue_mean_raw(1:length(sen_wue_mean_raw_agg),i) = sen_wue_mean_raw_agg;
            comp_sen_wue_agg_raw(1:length(sen_wue_agg_raw_agg),i) = sen_wue_agg_raw_agg;
            comp_std_wue_mean(1:length(std_wue_mean_agg),i) = std_wue_mean_agg;
            comp_std_wue_agg(1:length(std_wue_agg_agg),i) = std_wue_agg_agg;
            comp_sen_wue_up_mean(1:length(sen_wue_up_mean),i) = sen_wue_up_mean;
            comp_sen_wue_up_agg(1:length(sen_wue_up_agg),i) = sen_wue_up_agg;
            comp_sen_wue_low_mean(1:length(sen_wue_low_mean),i) = sen_wue_low_mean;
            comp_sen_wue_low_agg(1:length(sen_wue_low_agg),i) = sen_wue_low_agg;
            comp_sen_wue_up_mean_raw(1:length(sen_wue_up_mean_raw),i) = sen_wue_up_mean_raw;
            comp_sen_wue_up_agg_raw(1:length(sen_wue_up_agg_raw),i) = sen_wue_up_agg_raw;
            comp_sen_wue_low_mean_raw(1:length(sen_wue_low_mean_raw),i) = sen_wue_low_mean_raw;
            comp_sen_wue_low_agg_raw(1:length(sen_wue_low_agg_raw),i) = sen_wue_low_agg_raw;
            comp_pct_wue_mean(1:length(pct_wue_mean_agg),i) = pct_wue_mean_agg;
            comp_pct_wue_agg(1:length(pct_wue_agg_agg),i) = pct_wue_agg_agg;
        end
        comp_wue_mean(comp_wue_mean==0) = nan;
        comp_sen_wue_mean(comp_sen_wue_mean==0) = nan;
        comp_sen_wue_mean_raw(comp_sen_wue_mean_raw==0) = nan;
        comp_wue_agg(comp_wue_agg==0) = nan;
        comp_sen_wue_agg_raw(comp_sen_wue_agg_raw==0) = nan;
        comp_sen_wue_agg(comp_sen_wue_agg==0) = nan;
        comp_std_wue_mean(comp_std_wue_mean==0) = nan;
        comp_std_wue_agg(comp_std_wue_agg==0) = nan;
        comp_sen_wue_low_mean(comp_sen_wue_low_mean==0) = nan;
        comp_sen_wue_low_agg(comp_sen_wue_low_agg==0) = nan;
        comp_sen_wue_up_mean(comp_sen_wue_up_mean==0) = nan;
        comp_sen_wue_up_agg(comp_sen_wue_up_agg==0) = nan;
        comp_sen_wue_low_mean_raw(comp_sen_wue_low_mean_raw==0) = nan;
        comp_sen_wue_low_agg_raw(comp_sen_wue_low_agg_raw==0) = nan;
        comp_sen_wue_up_mean_raw(comp_sen_wue_up_mean_raw==0) = nan;
        comp_sen_wue_up_agg_raw(comp_sen_wue_up_agg_raw==0) = nan;
        comp_pct_wue_mean(comp_pct_wue_mean==0) = nan;
        comp_pct_wue_agg(comp_pct_wue_agg==0) = nan;
        clf;
        mcolor = [.5 1.0 .83; .90 .17 .31; .0 .0039 1.00; 1.00 .57 .69; .2 .5 0;...
            .00 .75 1.00; .76 .60 .42; .45 .63 .76; .53 .33 .04; 0 0 0];
        data_mean = comp_wue_mean;
        data_agg = comp_wue_agg;
        id = '(a)';
        yname = {wue_name; units};
        idx_plot = find(sum(isfinite(data_mean))~=0);
        idx = find(sum(isfinite(data_mean))==0);
        data_mean(:,idx) = [];
        idx = find(sum(isfinite(data_agg))==0);
        data_agg(:,idx) = [];
        idx = find(sum(isfinite(data_mean))==1);
        data_mean(:,idx) = nan;
        idx = find(sum(isfinite(data_agg))==1);
        data_agg(:,idx) = nan;
        [row, col] = size(data_agg);
        R2 = [];
        for j = 1:col
            data_agg(~isfinite(data_agg(:,j)),j) = nan;
            % R2 between uWUEhour and uWUE_season
            idx = isfinite(data_mean(:,j))&isfinite(data_agg(:,j));
            R2 = [R2; corr(data_mean(idx,j), data_agg(idx,j))^2];
        end
        subplot(2, 2, 1:2)
        h2 = violinplot(data_agg);
        hold on
        h1 = violinplot(data_mean);
        for j = 1:length(h1)
            h1(j).ShowData = 0;
            h1(j).ViolinColor = mcolor(2,:);
            h1(j).ShowMean = 0;
            h1(j).BoxColor = 'r';
            h1(j).BoxWidth = 0.01;
            h1(j).MedianPlot.Marker = '^';
            h1(j).MedianPlot.SizeData = 60;
            h1(j).MedianColor = 'none';
        end
        for j = 1:length(h2)
            h2(j).BoxColor = 'k';
            h2(j).BoxWidth = 0.04;
            h2(j).ViolinColor = 'k';
            h2(j).ShowData = 0;
            h2(j).EdgeColor = 'none';
            h2(j).ViolinAlpha = 0.2;
            h2(j).MedianPlot.SizeData = 20;
            h2(j).ShowMean = 0;
        end
        set(gca,'xticklabel',IGBP_name_w_yr(idx_plot))
        ylabel(yname, 'FontSize',12)
        t = text(0.02,1.1,[id],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top',...
            'FontSize',12);
        box on
        axis tight
        ylim(yrange)
        x_pos = linspace(0.02, 0.9, length(h1));
        for j = 1:length(h1)
            t = text(x_pos(j),0.95,['R^2=' num2str(sprintf('%1.2f',R2(j)))],...
                'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k');
        end
        ax1 = axes('Position',[0 0 1 1],'Visible','off');
        axes(ax1) % sets ax1 to current axes
        lgd = legend([h1(1).ViolinPlot, h2(1).ViolinPlot], {wue_med_lgd, wue_agg_lgd},...
            'orientation','horizontal', 'box','off');
        lgd.FontSize = 12;
        lgd.Position = [0.5 0.93 0.3643 0.0452];
        % scatter plot for normalized trend
        h_type = [];
        r2_type = [];
        mean_sen_range = [];
        agg_sen_range = [];
        wue_tmp = [];
        for i = 1:length(idx_plot)
            tmp_mean = comp_sen_wue_mean(:, idx_plot(i))*100;
            tmp_mean(~isfinite(tmp_mean)) = [];
            tmp_agg = comp_sen_wue_agg(:, idx_plot(i))*100;
            tmp_agg(~isfinite(tmp_agg)) = [];
            tmp_mean_low = comp_sen_wue_low_mean(:, idx_plot(i))*100;
            tmp_mean_low(~isfinite(tmp_mean_low)) = [];
            tmp_agg_low = comp_sen_wue_low_agg(:, idx_plot(i))*100;
            tmp_agg_low(~isfinite(tmp_agg_low)) = [];
            tmp_mean_up = comp_sen_wue_up_mean(:, idx_plot(i))*100;
            tmp_mean_up(~isfinite(tmp_mean_up)) = [];
            tmp_agg_up = comp_sen_wue_up_agg(:, idx_plot(i))*100;
            tmp_agg_up(~isfinite(tmp_agg_up)) = [];
            tmp = tmp_agg_up-tmp_agg_low;
            agg_sen_range = [agg_sen_range; tmp];
            tmp = tmp_mean_up-tmp_mean_low;
            mean_sen_range = [mean_sen_range; tmp];
            tmp = [tmp_mean, tmp_agg];
            wue_tmp = [wue_tmp; tmp];
            subplot(2, 2, 3)
            if (i==1)
                plot([-20 30], [-20 30], 'k--', 'linewidth', 1)
                ylabel([wue_agg_lgd ' trend (% yr^{-1})'], 'FontSize', 12)
                xlabel([wue_med_lgd ' trend (% yr^{-1})'], 'FontSize', 12)
            end
            hold on
            if (length(tmp_mean)>0 & length(tmp_agg)>0)
                hhh = errorbar(tmp_mean, tmp_agg, ...
                        tmp_agg_low, tmp_agg_up, tmp_mean_low, tmp_mean_up,...
                        'o', 'MarkerSize', 7,...
                        'Color', mcolor(i,:), 'MarkerFaceColor', mcolor(i,:));
                h_type = [h_type; hhh];
                r2_type = [r2_type; corr(tmp_mean, tmp_agg)^2];
            end
        end
        axis([-20 30 -20 30])
        % make an invisible axis handle
        ax_copy = axes('Position',get(gca,'Position'),'Visible','Off');
        tmp = corr(wue_tmp(:,1), wue_tmp(:,2))^2;
        t = text(0.02,0.98,['(b) R^2 = ' num2str(sprintf('%2.2f',tmp))],...
                'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                'FontSize',12);
        lgd = legend(ax_copy, [h_type],...
                    {char(IGBP_name(idx_plot))},...
                    'orientation','vertical', 'box','off');
        lgd.FontSize = 12;
        lgd.Position = [0.23 0.22 0.3643 0.0452];

        % density scatter for annual WUE pct changes
        subplot(2, 2, 4)
        wue_pct_mean = comp_pct_wue_mean(isfinite(comp_pct_wue_mean));
        wue_pct_agg = comp_pct_wue_agg(isfinite(comp_pct_wue_agg));
        ppp = dscatter(wue_pct_mean(:), wue_pct_agg(:));
        colormap('Hot')
        box on
        hold on
        tmp = [wue_pct_mean(:), wue_pct_agg(:)];
        tmp = tmp(isfinite(tmp));
        plot([min(tmp), max(tmp)], ...
            [min(tmp), max(tmp)], 'k--')
        % line fit
        rg_fit_tmp = polyfit(wue_pct_mean(:), wue_pct_agg(:), 1);
        xxx = linspace(min(tmp),max(tmp),10);
        fff = polyval(rg_fit_tmp,xxx);
        plot(xxx,fff,'b-','Linewidth',2)
        axis([min(tmp), max(tmp), min(tmp), max(tmp)])
        R2 = corr(wue_pct_mean(:), wue_pct_agg(:))^2;
        t = text(0.02,0.98,['(c) R^2 = ' num2str(sprintf('%1.2f',R2))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top',...
            'FontSize',12);
        ylabel([wue_agg_lgd ' IAV (% yr^{-1})'], 'FontSize', 12)
        xlabel([wue_med_lgd ' IAV (% yr^{-1})'], 'FontSize', 12)
        figname=['../plots_hourly/uWUE_trend_full_v7/WUE_dist_comp_' season_name '_' wue '.jpeg'];
        print('-dpng','-r300',[figname]);
    end
end

%% compare r values among mean WUE trends against trends inferred from each percentile
clc;
% analyze the top X% windows
pct_window = [1:99];
IGBP_included = unique(IGBP);% % exclude croplands and wetlands
idx = find(IGBP_included=='CRO'|IGBP_included=='WET');
IGBP_included(idx) = [];
plot_id = 1;
pct_best_r = nan(length(IGBP_included), 4, 4);
pct_best_r_site = nan(length(site_included), 4, 4);
site_count = 0;
for type = 1:4
    if (type==1) % grs, uwue
        wue_data = wue_hourly_agg;
        col_wue = 3;
        wue_name = '{\it uWUE}';
        units = '(g C hPa^{0.5} kg H_2O^{-1})';
        fig_header = 'uwue_pct_cor_grs_';
    elseif (type==2) % grs, iwue
        wue_data = wue_hourly_agg;
        col_wue = 15;
        wue_name = '{\it WUE_{ei}}';
        units = '(g C hPa kg H_2O^{-1})';
        fig_header = 'iwue_pct_cor_grs_';
    elseif (type==3) % grs, uwue
        wue_data = wue_hourly_summer;
        col_wue = 3;
        wue_name = '{\it uWUE}';
        units = '(g C hPa^{0.5} kg H_2O^{-1})';
        fig_header = 'uwue_pct_cor_sum_';
    elseif (type==4) % grs, iwue
        wue_data = wue_hourly_summer;
        col_wue = 15;
        wue_name = '{\it WUE_{ei}}';
        units = '(g C hPa kg H_2O^{-1})';
        fig_header = 'iwue_pct_cor_sum_';
    end
    for i = 1:length(IGBP_included)
        % for each IGBP
        idx_IGBP = find(IGBP==IGBP_included(i));
        % for sites under each IGBP
        clear site_list
        % declare the matrixs for trend analysis
        sen_agg = nan(length(idx_IGBP), 4, length(pct_window)+2);
        h_agg = nan(length(idx_IGBP), 4, length(pct_window)+2);
        norm_trend_agg = nan(length(idx_IGBP), 4, length(pct_window)+2);
        r_value = nan(length(idx_IGBP), length(pct_window));
        p_value = nan(length(idx_IGBP), length(pct_window));
        for variable = 1:1
            clf;
            annual_change = [];
            site_valid = [];
            for j = 1:length(idx_IGBP) % loop through sites
                idx = find(wue_data(:,4)==site_included(idx_IGBP(j)));
                data_tmp = wue_data(idx, :);
                % define growing season as GPP>0 & Tair>0
                idx = find(data_tmp(:, 5)>0 & data_tmp(:, 8)>0);
                grs_data = data_tmp(idx,:);
                grs_data(~isfinite(grs_data(:,col_wue)),:) = [];
                yr_range = unique(grs_data(:, 1));
                if (length(yr_range)>=10)
                    if (variable==1)
                        siteID = dir(site_included(idx_IGBP(j)), 5:10);
                        site_list{j} = {[siteID '(' int2str(length(yr_range)) ')']};
                    end
                    clear grs_sum MAH
                    % MAH calculation
                    for window = 1:length(pct_window)
                        % define MAH at the corresponding percentile
                        pct_value = prctile(grs_data(:, col_wue), pct_window(window));
                        for yr = 1:length(yr_range) % loop through years
                            idx = find(grs_data(:, 1)==yr_range(yr));
                            data_site_year = grs_data(idx, :);
                            idx = find(data_site_year(:, col_wue)>pct_value);
                            MAH(yr, window) = length(idx);
                        end
                    end
                    % growing season mean/sum
                    for yr = 1:length(yr_range) % loop through years
                        idx = find(grs_data(:, 1)==yr_range(yr));
                        data_site_year = grs_data(idx, :);
                        if (variable==2 | variable==4)
                            grs_sum(yr) = sum(data_site_year(:, col_wue),'omitnan');
                        elseif (variable==1 | variable==3)
                            grs_sum(yr) = median(data_site_year(:, col_wue),'omitnan');
                        end
                    end
                    % correlation coefficient against mean WUE
                    for window = 1:length(pct_window)
                        xtmp = MAH(:, window);
                        ytmp = grs_sum(:);
                        idx = isfinite(xtmp) & isfinite(ytmp);
                        if (sum(idx)>1)
                            if (window==1)
                                site_valid = [site_valid; j];
                            end
                            [r, p] = corrcoef(xtmp(idx), ytmp(idx));
                            r_value(j, window) = r(1,2);
                            p_value(j, window) = p(1,2); % if p<0.05, the corresponding r_value is significant
                        end
                    end
                    idx = find(r_value(j, :)==max(r_value(j, :)));
                    if (length(idx)>0)
                        pct_best_r_site(site_included(idx_IGBP(j)), variable, type) = pct_window(min(idx));
                    end
                end
            end
            if (sum(isfinite(r_value))>0)
                % specify the percentile resulting in the highest correlation on
                % average (Zscheischler et al 2016)
                if (length(site_valid)>1)
                    r_value_ave = mean(r_value(site_valid, :));
                    idx = find(r_value_ave==max(r_value_ave));
                    pct = pct_window(idx);
                else
                    pct = pct_site;
                end
                pct_best_r(i, variable, type) = pct;
            end
        end
    end
end

%% compare r values between seasonal uWUE and MAH inferred from each uWUE percentile
clc;
% analyze the top X% windows
pct_window = [1:99];
IGBP_included = unique(IGBP);% % exclude croplands
idx = find(IGBP_included=='CRO'|IGBP_included=='WET');
IGBP_included(idx) = [];
plot_id = 1;
pct_best_r_agg = nan(1, 4);
for type = 1:4
    if (type==1) % grs, uwue
        wue_data = wue_hourly_agg;
        col_wue = 3;
        wue_name = '{\it uWUE}';
        units = '(g C hPa^{0.5} kg H_2O^{-1})';
        fig_header = 'MAH_uwue_grs_';
        xname = '{\it {MAH(uWUE}}';
        wue_med_lgd = '{\it uWUE_{median}}';
    elseif (type==2) % grs, iwue
        wue_data = wue_hourly_agg;
        col_wue = 15;
        wue_name = '{\it WUE_{ei}}';
        units = '(g C hPa kg H_2O^{-1})';
        fig_header = 'MAH_iwue_grs_';
        xname = '{\it {MAH(WUE_{ei}}}';
        wue_med_lgd = '{\it WUE_{ei, median}}';
    elseif (type==3) % grs, uwue
        wue_data = wue_hourly_summer;
        col_wue = 3;
        wue_name = '{\it uWUE}';
        units = '(g C hPa^{0.5} kg H_2O^{-1})';
        fig_header = 'MAH_uwue_sum_';
        xname = '{\it {MAH(uWUE}}';
        wue_med_lgd = '{\it uWUE_{median}}';
    elseif (type==4) % grs, iwue
        wue_data = wue_hourly_summer;
        col_wue = 15;
        wue_name = '{\it WUE_{ei}}';
        units = '(g C hPa kg H_2O^{-1})';
        fig_header = 'MAH_iwue_sum_';
        xname = '{\it {MAH(WUE_{ei}}}';
        wue_med_lgd = '{\it WUE_{ei, median}}';
    end
    site_name_agg = [];
    sites = [];
    r_value_agg = [];
    p_value_agg = [];
    for i = 1:length(IGBP_included)
        % for each IGBP
        idx_IGBP = find(IGBP==IGBP_included(i));
        % for sites under each IGBP
        clear site_list
        % declare the matrixs for trend analysis
        sen_agg = nan(length(idx_IGBP), 4, length(pct_window)+2);
        h_agg = nan(length(idx_IGBP), 4, length(pct_window)+2);
        norm_trend_agg = nan(length(idx_IGBP), 4, length(pct_window)+2);
        r_value = nan(length(idx_IGBP), length(pct_window));
        p_value = nan(length(idx_IGBP), length(pct_window));
        for variable = 1:1
            clf;
            annual_change = [];
            site_valid = [];
            for j = 1:length(idx_IGBP) % loop through sites
                idx = find(wue_data(:,4)==site_included(idx_IGBP(j)));
                data_tmp = wue_data(idx, :);
                % define growing season as GPP>0 & Tair>0
                idx = find(data_tmp(:, 5)>0 & data_tmp(:, 8)>0);
                grs_data = data_tmp(idx,:);
                % remove nan WUE data points
                grs_data(~isfinite(grs_data(:,col_wue)),:) = [];
                if (sum(isfinite(grs_data(:,col_wue)))>0)
                    yr_range = unique(grs_data(:, 1));
                end
                if (variable==1)
                    siteID = dir(site_included(idx_IGBP(j)), 5:10);
                    site_list{j} = {[siteID '(' int2str(length(yr_range)) ')']};
                end
                if (length(yr_range)>=10)
                    site_valid = [site_valid; j];
                    clear grs_sum MAH
                    % MAH calculation
                    for window = 1:length(pct_window)
                        % define MAH at the corresponding percentile
                        pct_value = prctile(grs_data(:, col_wue), pct_window(window));
                        for yr = 1:length(yr_range) % loop through years
                            idx = find(grs_data(:, 1)==yr_range(yr));
                            data_site_year = grs_data(idx, :);
                            idx = find(data_site_year(:, col_wue)>pct_value);
                            MAH(yr, window) = length(idx);
                        end
                    end
                    % growing season mean/sum
                    for yr = 1:length(yr_range) % loop through years
                        idx = find(grs_data(:, 1)==yr_range(yr));
                        data_site_year = grs_data(idx, :);
                        if (variable==2 | variable==4)
                            grs_sum(yr) = sum(data_site_year(:, col_wue),'omitnan');
                        elseif (variable==1 | variable==3)
                            grs_sum(yr) = median(data_site_year(:, col_wue),'omitnan');
                        end
                    end
                    % correlation coefficient against mean WUE
                    for window = 1:length(pct_window)
                        xtmp = MAH(:, window);
                        ytmp = grs_sum(:);
                        idx = isfinite(xtmp) & isfinite(ytmp);
                        [r, p] = corrcoef(xtmp(idx), ytmp(idx));
                        r_value(j, window) = r(1,2);
                        p_value(j, window) = p(1,2); % if p<0.05, the corresponding r_value is significant
                    end
                    idx = find(r_value(j, :)==max(r_value(j, :)));
                    if (type==1)
                        pct_best_r_site(site_included(idx_IGBP(j)), variable) = pct_window(min(idx));
                    end
                end
            end
            if (length(site_valid))
                site_name_agg = [site_name_agg, site_list];
                r_value_agg = [r_value_agg; r_value];
                p_value_agg = [p_value_agg; p_value];
            end
        end
    end
    idx = find(isfinite(sum(r_value_agg, 2))==1);
    tmp = r_value_agg(idx, :);
    signal = p_value_agg(idx, :);
    site_info = site_name_agg(idx);
    r_value_ave = mean(r_value_agg(idx, :));
    idx = find(r_value_ave==max(r_value_ave));
    pct = pct_window(idx);
    % aggregated line plots
    clf;
    p = plot(pct_window, tmp, 'k-', 'linewidth', 1);
    for i = 1:length(p)
        p(i).Color(4) = 0.2;
    end
    hold on
    p_ens = shadedErrorBar(pct_window, mean(tmp), std(tmp),...
                    'lineprops', {'-','color',[5 128 0]/255}, 'transparent', 1);

    plot([pct pct], [-1 1], ':', 'color',[5 128 0]/255, 'linewidth', 2)
    axis([pct_window(1) pct_window(end) -1 1])
    ylabel(['Correlation between {\itMAH_{78-99}} and ' wue_med_lgd ' trends'], 'FontSize', 12)
    xlabel(['Most Actic Hours ({\itMAH}) defined at each ' wue_name ' percentile'], 'FontSize', 12)
    xticks([10:20:90])
    if (type==1|type==3)
        xticklabels({[xname '_{10^{th}})'], [xname '_{30^{th}})'], ...
            [xname '_{50^{th}})'], [xname '_{70^{th}})'], [xname '_{90^{th}})']})
    elseif (type==2|type==4)
        xticklabels({[xname '_{, 10^{th}})'], [xname '_{, 30^{th}})'], ...
            [xname '_{, 50^{th}})'], [xname '_{, 70^{th}})'], [xname '_{, 90^{th}})']})
    end
    lgd = legend([p_ens.mainLine, p(1)],...
                        {'Average', 'Site-specific'},...
                        'orientation','horizontal', 'box','off');
    lgd.FontSize = 12;
    lgd.Position = [0.55 0.93 0.3643 0.0452];
    if (type==1)
        t = text(0.02,0.08,['{\itMAH(uWUE_{x^{th}}, year) = \Sigma_{i=1}^{n} 1_{hourly uWUE(year)>uWUE_{x^{th}}}}'],'Units',...
                            'Normalized', 'VerticalAlignment', 'Top');
        text(78, -0.9, ['\it {MAH(uWUE_{78^{th}})}'], 'FontSize', 12, 'Color', [5 128 0]/255)
    elseif (type==2)
        t = text(0.02,0.08,['{\itMAH(WUE_{ei}_{, x^{th}}, year) = \Sigma_{i=1}^{n} 1_{hourly WUE_{ei}(year)>WUE_{ei}_{, x^{th}}}}'],'Units',...
                            'Normalized', 'VerticalAlignment', 'Top');
        text(pct, -0.9, ['\it {MAH(WUE_{ei}_{, ' int2str(pct) '^{th}})}'], 'FontSize', 12, 'Color', [5 128 0]/255)
    elseif (type==3)
        t = text(0.02,0.08,['{\itMAH(uWUE_{x^{th}}, year) = \Sigma_{i=1}^{n} 1_{hourly uWUE(year)>uWUE_{x^{th}}}}'],'Units',...
                            'Normalized', 'VerticalAlignment', 'Top');
        text(pct+3, -0.9, ['\it {MAH(uWUE_{' int2str(pct) '^{th}})}'], ...
            'FontSize', 12, 'Color', [5 128 0]/255, 'BackgroundColor','w')
    elseif (type==4)
        t = text(0.02,0.08,['{\itMAH(WUE_{ei}_{x^{th}}, year) = \Sigma_{i=1}^{n} 1_{hourly WUE_{ei}(year)>WUE_{ei}_{x^{th}}}}'],'Units',...
                            'Normalized', 'VerticalAlignment', 'Top');
        text(pct+3, -0.9, ['\it {MAH(WUE_{ei}_{, ' int2str(pct) '^{th}})}'], ...
            'FontSize', 12, 'Color', [5 128 0]/255, 'BackgroundColor','w')
    end
    figname=['../plots_hourly/uWUE_trend_full_v7/' ...
                fig_header '_line_v2.jpeg'];
    print('-dpng','-r300',figname);
end

%% compare the WUE(ave) trends with and without MAHs. include CMIP6 results
id_name = ["(a)", "(b)", "(c)","(d)", "(e)", "(f)", "(g)", "(h)", "(i)", ...
    "(j)", "(k)", "(l)","(m)", "(n)", "(o)", "(p)", "(q)", "(r)", "(s)", "(t)", "(u)"];
clc;
% analyze the top X% windows
pct_window = [95];
IGBP_included = unique(IGBP);% % exclude croplands
idx = find(IGBP_included=='CRO'|IGBP_included=='WET');
IGBP_included(idx) = [];
for type = 1:4
    clear r_agg wue_cmip6
    CMIP6_LAT = [];
    CMIP6_LON = [];
    if (type==1) % grs, uwue
        wue_data = wue_hourly_agg;
        col_wue = 3;
        wue_name = '{\it uWUE}';
        units = '(g C hPa^{0.5} kg H_2O^{-1})';
        fig_header = 'uWUE_MAH_trend_grs';
        xname = '{\it {MAH(uWUE}}';
        wue_med_lgd = '{\it uWUE_{median}}';
        wue_agg_lgd = '{\it uWUE_{aggregate}}';
    elseif (type==2) % grs, iwue
        wue_data = wue_hourly_agg;
        col_wue = 15;
        wue_name = '{\it WUE_{ei}}';
        units = '(g C hPa kg H_2O^{-1})';
        fig_header = 'iWUE_MAH_trend_grs';
        xname = '{\it {MAH(WUE_{ei}}}';
        wue_med_lgd = '{\it WUE_{ei, median}}';
    elseif (type==3) % grs, uwue
        wue_data = wue_hourly_summer;
        col_wue = 3;
        wue_name = '{\it uWUE}';
        units = '(g C hPa^{0.5} kg H_2O^{-1})';
        fig_header = 'uWUE_MAH_trend_sum';
        xname = '{\it {MAH(uWUE}}';
        wue_med_lgd = '{\it uWUE_{median}}';
    elseif (type==4) % grs, iwue
        wue_data = wue_hourly_summer;
        col_wue = 15;
        wue_name = '{\it WUE_{ei}}';
        units = '(g C hPa kg H_2O^{-1})';
        fig_header = 'iWUE_MAH_trend_sum';
        xname = '{\it {MAH(WUE_{ei}}}';
        wue_med_lgd = '{\it WUE_{ei, median}}';
    end
    r_Sen_MAH_pct_type = [];
    r_Sen_MAH_pct_agg = [];
    wue_norm_trend = [];
    wue_agg_norm_trend = [];
    wue_norm_trend_type = [];
    wue_norm_trend_agg = [];
    wue_ts = [];
    wue_ts_type = [];
    wue_ts_agg = [];
    wue_norm_trend_p_MAH = [];
    wue_norm_trend_n_MAH = [];
    wue_norm_trend_agg_p_MAH = [];
    wue_norm_trend_agg_n_MAH = [];
    wue_agg_norm_trend_agg_p_MAH = [];
    wue_agg_norm_trend_agg_n_MAH = [];
    % read in CMIP6 WUE
    files = folderFiles(['../CMIP6_uWUE/'],'uwue*');
    [row, col] = size(files);
    for i = 1:row
        filename = strtrim(['../CMIP6_uWUE/' files(i,:)]);
        T = importdata(filename,',');
        wue_cmip6(:,:,i) = T;
    end
    % exclude ZA-Kru
    wue_cmip6(:,38,:) = [];
    tmp = size(wue_cmip6);
    cmip6_model = tmp(3);
    cmip6_yr = 1990:2014;
    cmip6_wue_trend_agg = [];
    cmip6_wue_trend_agg_p_MAH = [];
    cmip6_wue_trend_agg_n_MAH = [];
    cmip6_wue_trend_agg_w_MAH = [];
    count = 0;
    for i = 1:length(IGBP_included)
        plot_id = [];  
        % for each IGBP
        idx_IGBP = find(IGBP==IGBP_included(i));
        % for sites under each IGBP
        clear site_list
        % declare the matrixs for trend analysis
        sen_agg = nan(length(idx_IGBP), 4, length(pct_window)+2);
        h_agg = nan(length(idx_IGBP), 4, length(pct_window)+2);
        norm_trend_agg = nan(length(idx_IGBP), 4, length(pct_window)+2);
        MAH_trend_agg = nan(length(idx_IGBP), 4);
        for variable = 1:1
            clf;
            pct_threshold = pct_best_r(i, variable, type);
            pct_threshold_agg = pct_best_r_agg(variable, type);
            annual_change = [];
            for j = 1:length(idx_IGBP) % loop through sites
                idx = find(wue_data(:,4)==site_included(idx_IGBP(j)));
                data_tmp = wue_data(idx, :);
                % define growing season as GPP>0 & Tair>0
                idx = find(data_tmp(:, 5)>0 & data_tmp(:, 8)>0);
                grs_data = data_tmp(idx,:);
                % find the pct values
                pct_value = prctile(grs_data(:,col_wue), pct_threshold);
                pct_value_agg = prctile(grs_data(:,col_wue), pct_threshold_agg);
                yr_range = unique(grs_data(:, 1));
                if (length(yr_range)>=10)
                    plot_id = [plot_id; j];
                    idx_cmip6_month = [];
                    clear grs_sum MAH MAH_agg wue_detrend_MAH_type wue_detrend_MAH_agg wue_raw wue_agg
                    for yr = 1:length(yr_range) % loop through years
                        idx = find(grs_data(:, 1)==yr_range(yr));
                        data_site_year = grs_data(idx, :);
                        tmp = find(mon_sum>=min(data_site_year(:,2))&...
                            mon_sum<=max(data_site_year(:,2)));
                        year = find(yr_range(yr)==cmip6_yr);
                        idx_cmip6_month = [idx_cmip6_month; tmp'+(year-1)*12];
                        tmp = data_site_year(:, col_wue);
                        wue_raw(yr) = median(tmp,'omitnan');
                        wue_agg(yr) = sum(data_site_year(:,5),'omitnan')*...
                            (mean(data_site_year(:,7),'omitnan')^.5)/...
                            sum(data_site_year(:,6),'omitnan');
                        idx = find(data_site_year(:, col_wue)>=pct_value);
                        tmp(idx) = nan;
                        wue_detrend_MAH_type(yr) = median(tmp,'omitnan');
                        MAH(yr) = length(idx);
                        tmp = data_site_year(:, col_wue);
                        idx = find(data_site_year(:, col_wue)>=pct_value_agg);
                        tmp(idx) = nan;
                        wue_detrend_MAH_agg(yr) = median(tmp,'omitnan');
                        MAH_agg(yr) = length(idx);
                    end
                        if (variable==1)
                            siteID = dir(site_included(idx_IGBP(j)), 5:10);
                            site_list{j} = {[siteID '(' int2str(length(yr_range)) ')']};
                        end
                        tmp = wue_raw;
                        idx_valid = find(isfinite(tmp)==1);
                        if (length(idx_valid)>=10)
                            CMIP6_LAT = [CMIP6_LAT; LAT((idx_IGBP(j)))];
                            CMIP6_LON = [CMIP6_LON; LONG((idx_IGBP(j)))];
                            count = count+1;
                            % cmip6 wue at the corresponding years
                            cmip6_wue = [];
                            cmip6_wue_trend = [];
                            for model = 1:cmip6_model
                                % match the measurement period
                                xxx = squeeze(wue_cmip6(:,count, model));
                                yyy = nan(length(xxx),1);
                                yyy(idx_cmip6_month) = xxx(idx_cmip6_month); 
                                tmp = mean(reshape(yyy, 12, length(cmip6_yr)), 'omitnan');
                                if (model==1)
                                    idx = [];
                                    for yr = 1:length(yr_range)
                                        idx = [idx; find(cmip6_yr==yr_range(yr))];
                                    end
                                end
                                model_wue = (tmp(idx)-mean(tmp(idx),'omitnan'))/mean(tmp(idx),'omitnan');
                                cmip6_wue = [cmip6_wue, model_wue'];
                                model_trend = diff(tmp(idx))./diff(yr_range')/mean(tmp(idx),'omitnan');
                                cmip6_wue_trend = [cmip6_wue_trend, model_trend'];
                            end
                            cmip6_wue_trend_nt(count) = length(yr_range)-1;
                            cmip6_wue_trend_agg = [cmip6_wue_trend_agg; cmip6_wue_trend];
                            datain = nan(length(idx_valid),4);
                            datain(:,1) = yr_range(idx_valid);
                            datain(:,2) = (wue_raw(idx_valid))/mean(wue_raw(idx_valid));
                            datain(:,3) = (wue_detrend_MAH_type(idx_valid)-mean(wue_detrend_MAH_type(idx_valid)))...
                                /mean(wue_detrend_MAH_type(idx_valid));
                            datain(:,4) = (wue_detrend_MAH_agg(idx_valid)-mean(wue_detrend_MAH_agg(idx_valid)))...
                                /mean(wue_detrend_MAH_agg(idx_valid));
                            datain(:,5) = MAH_agg(idx_valid);
                            datain(:,6) = (wue_agg(idx_valid))/mean(wue_agg(idx_valid));
                            wue_ts = [wue_ts; diff(datain(:,2))./diff(datain(:,1))/mean(datain(:,2))];
                            wue_ts_type = [wue_ts_type; diff(datain(:,3))./diff(datain(:,1))/mean(datain(:,3))];
                            wue_ts_agg = [wue_ts_agg; diff(datain(:,4))./diff(datain(:,1))/mean(datain(:,4))];
                            % MAH trends
                            [taub tau h sig Z S sigma sen_MAH n senplot CIlower CIupper D Dall C3 nsigma] =...
                                ktaub(datain(:, [1, 5]), 0.05);
                            % WUE trends, CMIP6
                            cmip6_tmp = [];
                            for model = 1:cmip6_model
                                tmp = cmip6_wue(idx_valid,model);
                                if (sum(isfinite(tmp))>2)
                                    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] =...
                                        ktaub([datain(:,1), tmp], 0.05);
                                else
                                    sen = NaN;
                                end
                                cmip6_tmp = [cmip6_tmp, sen];
                            end
                            cmip6_wue_trend_agg_w_MAH = [cmip6_wue_trend_agg_w_MAH; cmip6_tmp];
                            % WUE trends, observations
                            [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] =...
                                ktaub(datain(:, [1, 2]), 0.05);
                            wue_norm_trend = [wue_norm_trend; sen];   
                            % WUE trends, agg, observations
                            [taub tau h sig Z S sigma sen_agg n senplot CIlower CIupper D Dall C3 nsigma] =...
                                ktaub(datain(:, [1, 6]), 0.05);
                            wue_agg_norm_trend = [wue_agg_norm_trend; sen_agg];   
                            if (sen_MAH>0)
                                wue_norm_trend_p_MAH = [wue_norm_trend_p_MAH; sen];
                                wue_agg_norm_trend_agg_p_MAH = [wue_agg_norm_trend_agg_p_MAH; sen_agg];
                                cmip6_tmp = [];
                                for model = 1:cmip6_model
                                    tmp = cmip6_wue(idx_valid,model);
                                    if (sum(isfinite(tmp))>2)
                                        [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] =...
                                            ktaub([datain(:,1), tmp], 0.05);
                                    else
                                        sen = NaN;
                                    end
                                    cmip6_tmp = [cmip6_tmp, sen];
                                end
                                cmip6_wue_trend_agg_p_MAH = [cmip6_wue_trend_agg_p_MAH; cmip6_tmp];
                            elseif (sen_MAH<0)
                                wue_norm_trend_n_MAH = [wue_norm_trend_n_MAH; sen];
                                wue_agg_norm_trend_agg_n_MAH = [wue_agg_norm_trend_agg_n_MAH; sen_agg];
                                cmip6_tmp = [];
                                for model = 1:cmip6_model
                                    tmp = cmip6_wue(idx_valid,model);
                                    if (sum(isfinite(tmp))>2)
                                        [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] =...
                                            ktaub([datain(:,1), tmp], 0.05);
                                    else
                                        sen = NaN;
                                    end
                                    cmip6_tmp = [cmip6_tmp, sen];
                                end
                                cmip6_wue_trend_agg_n_MAH = [cmip6_wue_trend_agg_n_MAH; cmip6_tmp];
                            end
                            % WUE trends without WUE during type-specific MAHs
                            [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] =...
                                ktaub(datain(:, [1, 3]), 0.05);
                            wue_norm_trend_type = [wue_norm_trend_type; sen];
                            % WUE trends without WUE during MAHs
                            [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] =...
                                ktaub(datain(:, [1, 4]), 0.05);
                            wue_norm_trend_agg = [wue_norm_trend_agg; sen];
                            if (sen_MAH>0)
                                wue_norm_trend_agg_p_MAH = [wue_norm_trend_agg_p_MAH; sen];
                            elseif (sen_MAH<0)
                                wue_norm_trend_agg_n_MAH = [wue_norm_trend_agg_n_MAH; sen];
                            end
                        end
                    end                
                end
            end
    end
    % compare WUE trend with and without MAH
    close all;
    wue_ts = wue_ts*100;
    wue_ts_agg = wue_ts_agg*100;
    subplot(3, 1, 1)
    tmp = [wue_norm_trend, wue_agg_norm_trend, ...
        mean(cmip6_wue_trend_agg_w_MAH,2,'omitnan')]*100;
    hhh = violinplot(tmp);
    median_value = median(tmp);
    std_value = std(tmp);
    box on
    for j = 1:length(hhh)
        hhh(j).BoxColor = 'k';
        hhh(j).BoxWidth = 0.04;
        hhh(j).ShowData = 0;
        hhh(j).ViolinColor = mcolor(j+1,:);
        if (j==length(hhh))
            hhh(j).ViolinColor = 'k';
        end
        hhh(j).ShowData = 1;
    end
    hold on
    plot([0.5 length(hhh)+.5], [0, 0], '--', 'Color', 'k')
    ylim([-4 10])
    t = text(0.02,1.2,['(a) When {\it MAH_{uWUE^{78^{th}-99^{th}}}} trends are positive or negative (' int2str(length(wue_norm_trend)) ' sites)'],...
        'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k', 'FontSize', 12);
    x_pos = linspace(0.12, 0.8, length(hhh));
    for j = 1:length(hhh)
        t = text(x_pos(j),0.95,[num2str(sprintf('%2.1f',median_value(j))) char(177) ...
            num2str(sprintf('%2.1f',std_value(j)))],...
            'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k');
    end
    set(gca,'xticklabel',...
        {['Observed ' wue_med_lgd], ['Observed ' wue_agg_lgd], 'CMIP6'},...
        'FontSize', 12)
    ylabel([wue_name ' trends (% yr^{-1})'], 'FontSize', 12)
    subplot(3, 1, 2)
    tmp = [wue_norm_trend_p_MAH, wue_agg_norm_trend_agg_p_MAH, ...
        mean(cmip6_wue_trend_agg_p_MAH,2,'omitnan')]*100;
    hhh = violinplot(tmp);
    median_value = median(tmp);
    std_value = std(tmp);
    box on
    for j = 1:length(hhh)
        hhh(j).BoxColor = 'k';
        hhh(j).BoxWidth = 0.04;
        hhh(j).ShowData = 0;
        hhh(j).ViolinColor = mcolor(j+1,:);
        if (j==length(hhh))
            hhh(j).ViolinColor = 'k';
        end
        hhh(j).ShowData = 1;
    end
    hold on
    plot([0.5 length(hhh)+.5], [0, 0], '--', 'Color', 'k')
    ylim([-4 10])
    t = text(0.02,1.2,['(b) When {\it MAH_{uWUE^{78^{th}-99^{th}}}} trends are positive (' int2str(length(wue_norm_trend_p_MAH)) ' sites)'],...
        'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k', 'FontSize', 12);
    for j = 1:length(hhh)
        t = text(x_pos(j),0.95,[num2str(sprintf('%2.1f',median_value(j))) char(177) ...
            num2str(sprintf('%2.1f',std_value(j)))],...
            'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k');
    end
    set(gca,'xticklabel',...
        {['Observed ' wue_med_lgd], ['Observed ' wue_agg_lgd], 'CMIP6'},...
        'FontSize', 12)
    ylabel([wue_name ' trends (% yr^{-1})'], 'FontSize', 12)
    subplot(3, 1, 3)
    tmp = [wue_norm_trend_n_MAH, wue_agg_norm_trend_agg_n_MAH, ...
        mean(cmip6_wue_trend_agg_n_MAH,2,'omitnan')]*100;
    hhh = violinplot(tmp);
    median_value = median(tmp);
    std_value = std(tmp);
    box on
    for j = 1:length(hhh)
        hhh(j).BoxColor = 'k';
        hhh(j).BoxWidth = 0.04;
        hhh(j).ShowData = 0;
        hhh(j).ViolinColor = mcolor(j+1,:);
        if (j==length(hhh))
            hhh(j).ViolinColor = 'k';
        end
        hhh(j).ShowData = 1;
    end
    hold on
    plot([0.5 length(hhh)+.5], [0, 0], '--', 'Color', 'k')
    ylim([-4 10])
    t = text(0.02,1.2,['(c) When {\it MAH_{uWUE^{78^{th}-99^{th}}}} trends are negative (' int2str(length(wue_norm_trend_n_MAH)) ' sites)'],...
        'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k', 'FontSize', 12);
    for j = 1:length(hhh)
        t = text(x_pos(j),0.95,[num2str(sprintf('%2.1f',median_value(j))) char(177) ...
            num2str(sprintf('%2.1f',std_value(j)))],...
            'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k');
    end
    set(gca,'xticklabel',...
        {['Observed ' wue_med_lgd], ['Observed ' wue_agg_lgd], 'CMIP6'},...
        'FontSize', 12)
    ylabel([wue_name ' trends (% yr^{-1})'], 'FontSize', 12)

    figname=['../plots_hourly/uWUE_trend_full_v7/' ...
                        fig_header '.jpeg'];
    print('-dpng','-r300',figname);
end

%% MAH sensitivity to WUE threshold; MAH correlation to WUE
% compare the WUE trends against # of data points above type specific WUE pct
id_name = ["(a)", "(b)", "(c)","(d)", "(e)", "(f)", "(g)", "(h)", "(i)", ...
    "(j)", "(k)", "(l)","(m)", "(n)", "(o)", "(p)", "(q)", "(r)", "(s)", "(t)", "(u)"];
clc;
close all;
for type = 1:4
    if (type==1) % grs, uwue
        wue_data = wue_hourly_agg;
        col_wue = 3;
        wue_name = '{\it uWUE}';
        units = '(gC hPa^{0.5} kg H_2O^{-1})';
        fig_header = 'cor_uWUE_MAH_misc_grs';
        xname = '{\it {MAH(uWUE}}';
        wue_med_lgd = '{\it uWUE_{median}}';
        wue_agg_lgd = '{\it uWUE_{aggregate}}';
    elseif (type==2) % grs, iwue
        wue_data = wue_hourly_agg;
        col_wue = 15;
        wue_name = '{\it WUE_{ei}}';
        units = '(gC hPa kg H_2O^{-1})';
        fig_header = 'cor_iWUE_MAH_misc_grs';
        xname = '{\it {MAH(WUE_{ei}}}';
        wue_med_lgd = '{\it WUE_{ei, median}}';
    elseif (type==3) % grs, uwue
        wue_data = wue_hourly_summer;
        col_wue = 3;
        wue_name = '{\it uWUE}';
        units = '(gC hPa^{0.5} kg H_2O^{-1})';
        fig_header = 'cor_uWUE_MAH_misc_sum';
        xname = '{\it {MAH(uWUE}}';
        wue_med_lgd = '{\it uWUE_{median}}';
    elseif (type==4) % grs, iwue
        wue_data = wue_hourly_summer;
        col_wue = 15;
        wue_name = '{\it WUE_{ei}}';
        units = '(gC hPa kg H_2O^{-1})';
        fig_header = 'cor_iWUE_MAH_misc_sum';
        xname = '{\it {MAH(WUE_{ei}}}';
        wue_med_lgd = '{\it WUE_{ei, median}}';
    end
    clear r_agg
    % analyze the top X% windows
    pct_window = [95];
    IGBP_included = unique(IGBP);% % exclude croplands
    idx = find(IGBP_included=='CRO'|IGBP_included=='WET');
    IGBP_included(idx) = [];
    r_Sen_MAH_pct_type = [];
    r_Sen_MAH_pct_agg = [];
    r_Sen_MAH_pct_site = [];
    r_Sen_MAH_pct_base = [];
    r_Sen_MAH_pct_base_v2 = [];
    r_sup_agg = [];
    site_MAH_ana = [];
    for i = 1:length(IGBP_included)
        plot_id = [];
        if (IGBP_included(i)=='CRO')
            class = 'Croplands';
            pct_threshold = 100;
        elseif (IGBP_included(i)=='CSH')
            class = 'Closed Shrublands';
            pct_threshold = 100;
        elseif (IGBP_included(i)=='DBF')
            class = 'Deciduous Broadleaf Forests';
            pct_threshold = 88;
        elseif (IGBP_included(i)=='EBF')
            class = 'Evergreen Broadleaf Forests';
            pct_threshold = 93;
        elseif (IGBP_included(i)=='ENF')
            class = 'Evergreen Needleleaf Forests';
            pct_threshold = 76;
        elseif (IGBP_included(i)=='GRA')
            class = 'Grasslands';
            pct_threshold = 69;
        elseif (IGBP_included(i)=='MF')
            class = 'Mixed Forests';
            pct_threshold = 87;
        elseif (IGBP_included(i)=='OSH')
            class = 'Open Shrublands';
            pct_threshold = 100;
        elseif (IGBP_included(i)=='SAV')
            class = 'Savannas';
            pct_threshold = 100;
        elseif (IGBP_included(i)=='WET')
            class = 'Permanent Wetlands';
            pct_threshold = 100;
        elseif (IGBP_included(i)=='WSA')
            class = 'Woody Savannas';
            pct_threshold = 84;
        end
        % for each IGBP
        idx_IGBP = find(IGBP==IGBP_included(i));
        % for sites under each IGBP
        clear site_list
        % declare the matrixs for trend analysis
        sen_agg = nan(length(idx_IGBP), 4, length(pct_window)+2);
        h_agg = nan(length(idx_IGBP), 4, length(pct_window)+2);
        norm_trend_agg = nan(length(idx_IGBP), 4, length(pct_window)+2);
        MAH_trend_agg = nan(length(idx_IGBP), 4);
        for variable = 1:1
            annual_change = [];
            count = 0;
            for j = 1:length(idx_IGBP) % loop through sites
                if (site_included(idx_IGBP(j))<166)
                    pct_threshold_type = pct_best_r(i, variable, type);
                    pct_threshold_site = pct_best_r_site(site_included(idx_IGBP(j)), variable, type);
                    pct_threshold_agg = pct_best_r_agg(type);
                    idx = find(wue_data(:,4)==site_included(idx_IGBP(j)));
                    data_tmp = wue_data(idx, :);
                    % define growing season as GPP>0 & Tair>0
                    idx = find(data_tmp(:, 5)>0 & data_tmp(:, 8)>0);
                    grs_data = data_tmp(idx,:);
                    % find the pct values
                    pct_value_type = prctile(grs_data(:,col_wue), pct_threshold_type);
                    pct_value_site = prctile(grs_data(:,col_wue), pct_threshold_site);
                    pct_value_agg = prctile(grs_data(:,col_wue), pct_threshold_agg);
                    pct_value_base = prctile(grs_data(:,col_wue), 50);
                    pct_value_base_low = prctile(grs_data(:,col_wue), 40);
                    pct_value_base_up = prctile(grs_data(:,col_wue), 60);
                    yr_range = unique(grs_data(:, 1));
                    if (length(yr_range)>=10)
                        plot_id = [plot_id; j];
                        count = count+1;
                        clear grs_sum MAH_site MAH_type MAH_agg MAH_base MAH_base_v2
                        sup = [];
                        for yr = 1:length(yr_range) % loop through years
                            idx = find(grs_data(:, 1)==yr_range(yr));
                            data_site_year = grs_data(idx, :);
                            for window = 1:length(pct_window)
                                % use different color to indicate different trends
                                % calculated by different metrics
                                tmp = prctile(data_site_year(:, col_wue),pct_window(window));
                                if (pct_window(window)>=50)
                                    idx = find(data_site_year(:, col_wue)>=tmp);
                                else
                                    idx = find(data_site_year(:, col_wue)<=tmp);
                                end
                                if (variable==2 | variable==4)
                                    grs_sum(yr, window) = sum(data_site_year(idx, col_wue));
                                elseif (variable==1 | variable==3)
                                    grs_sum(yr, window) = median(data_site_year(idx, col_wue));
                                end
                            end
                            if (variable==2 | variable==4)
                                grs_sum(yr, length(pct_window)+1) = sum(data_site_year(:, col_wue),'omitnan');
                            elseif (variable==1 | variable==3)
                                grs_sum(yr, length(pct_window)+1) = median(data_site_year(:, col_wue),'omitnan');
                            end
                            tmp = [sum(data_site_year(:, 5),'omitnan'), sum(data_site_year(:, 6),'omitnan'),...
                                mean(data_site_year(:, 7),'omitnan'), mean(data_site_year(:, 9),'omitnan')];
                            sup = [sup; tmp];
                            idx = find(data_site_year(:, col_wue)>=pct_value_type);
                            MAH_type(yr) = length(idx);
                            idx = find(data_site_year(:, col_wue)>=pct_value_site);
                            MAH_site(yr) = length(idx);
                            idx = find(data_site_year(:, col_wue)>=pct_value_agg);
                            MAH_agg(yr) = length(idx);
                            idx = find(data_site_year(:, col_wue)>=pct_value_base_low & data_site_year(:, col_wue)<=pct_value_base_up);
                            MAH_base(yr) = length(idx);
                            idx = find(data_site_year(:, col_wue)>=pct_value_base & data_site_year(:, col_wue)<pct_value_agg);
                            MAH_base_v2(yr) = length(idx);
                        end
                        annual_change = [annual_change; diff(grs_sum)./mean(grs_sum,'omitnan')];
                    if (variable==1)
                        siteID = dir(site_included(idx_IGBP(j)), 5:10);
                        site_list{j} = {[siteID '(' int2str(length(yr_range)) ')']};
                    end
                        % plot trends and MAH distributions for each site
                        tmp = grs_sum(:, length(pct_window)+1);
                        idx = find(isfinite(tmp)==1);
                        if (length(idx)>=10)
                            datain = nan(length(idx),2);
                            datain(:,1) = yr_range(idx);
                            datain(:,2) = (tmp(idx)-mean(tmp(idx)))/mean(tmp(idx));
                            r = corr(datain(:,2), MAH_type(idx)');
                            r_Sen_MAH_pct_type = [r_Sen_MAH_pct_type; r];
                            r = corr(datain(:,2), MAH_site(idx)');
                            r_Sen_MAH_pct_site = [r_Sen_MAH_pct_site; r];
                            r = corr(datain(:,2), MAH_agg(idx)');
                            r_Sen_MAH_pct_agg = [r_Sen_MAH_pct_agg; r];
                            r = corr(datain(:,2), MAH_base(idx)');
                            r_Sen_MAH_pct_base = [r_Sen_MAH_pct_base; r];
                            r = corr(datain(:,2), MAH_base_v2(idx)');
                            r_Sen_MAH_pct_base_v2 = [r_Sen_MAH_pct_base_v2; r];
                            sup_tmp = sup(idx,:);
                            tmp = [];
                            for kk = 1:4
                                idx_val = isfinite(datain(:,2))&isfinite(sup_tmp(:,kk));
                                tmp = [tmp, corr(datain(idx_val,2), sup_tmp(idx_val,kk))];
                            end
                            r_sup_agg = [r_sup_agg; tmp];
                            site_MAH_ana = [site_MAH_ana; site_included(idx_IGBP(j))];
                            [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] =...
                                ktaub(datain, 0.05);
                        end
                    end
                end
            end
        end
    end
    % violin plot for corr distribution between MAH and WUE trend
    clf;
    tmp = [r_Sen_MAH_pct_site, r_Sen_MAH_pct_type, r_Sen_MAH_pct_agg, ...
        r_Sen_MAH_pct_base, r_Sen_MAH_pct_base_v2, r_sup_agg];
    median_value = median(tmp, 'omitnan');
    std_value = std(tmp, 'omitnan');
    hhh = violinplot(tmp);
    box on
    for j = 1:length(hhh)
        hhh(j).BoxColor = 'k';
        hhh(j).BoxWidth = 0.04;
        hhh(j).ShowData = 0;
        if (j<=5)
            hhh(j).ViolinColor = mcolor(2,:);
        else
            hhh(j).ViolinColor = mcolor(3,:);
        end
        hhh(j).ShowData = 1;
    end
    x_pos = linspace(0.03, 0.88, length(hhh));
    for j = 1:length(hhh)
        t = text(x_pos(j),0.05,[num2str(sprintf('%2.1f',median_value(j))) char(177) ...
            num2str(sprintf('%2.1f',std_value(j)))],...
            'Units','Normalized', 'VerticalAlignment', 'Top', 'color', 'k');
    end
    if (type==1)
        set(gca,'xticklabel',...
            {'{\it MAH_{site}}', '{\it MAH_{eco}}', ...
            '{\it MAH_{uWUE^{78^{th}-99^{th}}}}', '{\it MAH_{uWUE^{40^{th}-60^{th}}}}', ...
            '{\it MAH_{uWUE^{50^{th}-77^{th}}}}',...
            'GPP', 'ET', 'VPD', 'CO_2'},...
            'FontSize', 12, 'XTickLabelRotation', -35)
    elseif (type==2)
        set(gca,'xticklabel',...
            {'{\it MAH_{site}}', '{\it MAH_{eco}}', ...
            '{\it MAH_{WUE_{ei}^{78^{th}-99^{th}}}}', '{\it MAH_{WUE_{ei}^{40^{th}-60^{th}}}}', ...
            '{\it MAH_{WUE_{ei}^{50^{th}-77^{th}}}}',...
            'GPP', 'ET', 'VPD', 'CO_2'},...
            'FontSize', 12, 'XTickLabelRotation', -35)
    elseif (type==3)
        set(gca,'xticklabel',...
            {'{\it MAH_{site}}', '{\it MAH_{eco}}', ...
            '{\it MAH_{uWUE^{80^{th}-99^{th}}}}', '{\it MAH_{uWUE^{40^{th}-60^{th}}}}', ...
            '{\it MAH_{uWUE^{50^{th}-79^{th}}}}',...
            'GPP', 'ET', 'VPD', 'CO_2'},...
            'FontSize', 12, 'XTickLabelRotation', -35)
    elseif (type==4)
        set(gca,'xticklabel',...
            {'{\it MAH_{site}}', '{\it MAH_{eco}}', ...
            '{\it MAH_{WUE_{ei}^{83^{rd}-99^{th}}}}', '{\it MAH_{WUE_{ei}^{40^{th}-60^{th}}}}', ...
            '{\it MAH_{WUE_{ei}^{50^{th}-82^{nd}}}}',...
            'GPP', 'ET', 'VPD', 'CO_2'},...
            'FontSize', 12, 'XTickLabelRotation', -35)
    end
    ylabel(['Correlation with ' wue_med_lgd ' trend'], 'FontSize', 12)
    figname=['../plots_hourly/uWUE_trend_full_v7/' ...
                fig_header '.jpeg'];
    print('-dpng','-r300',figname);
end

%% MAH sensitivity to WUE magnitude and PFT
id_name = ["(a)", "(b)", "(c)","(d)", "(e)", "(f)", "(g)", "(h)", "(i)", ...
    "(j)", "(k)", "(l)","(m)", "(n)", "(o)", "(p)", "(q)", "(r)", "(s)", "(t)", "(u)"];
clc;
close all
for type = 1:4
    if (type==1) % grs, uwue
        wue_data = wue_hourly_agg;
        col_wue = 3;
        wue_name = '{\it uWUE}';
        units = '(g C hPa^{0.5} kg H_2O^{-1})';
        fig_header = 'uwue_grs_MAH';
        xname = '{\it {MAH(uWUE}}';
        wue_med_lgd = '{\it uWUE_{median}}';
        wue_agg_lgd = '{\it uWUE_{aggregate}}';
    elseif (type==2) % grs, iwue
        wue_data = wue_hourly_agg;
        col_wue = 15;
        wue_name = '{\it WUE_{ei}}';
        units = '(g C hPa kg H_2O^{-1})';
        fig_header = 'iwue_grs_MAH';
        xname = '{\it {MAH(WUE_{ei}}}';
        wue_med_lgd = '{\it WUE_{ei, median}}';
    elseif (type==3) % grs, uwue
        wue_data = wue_hourly_summer;
        col_wue = 3;
        wue_name = '{\it uWUE}';
        units = '(g C hPa^{0.5} kg H_2O^{-1})';
        fig_header = 'uwue_sum_MAH';
        xname = '{\it {MAH(uWUE}}';
        wue_med_lgd = '{\it uWUE_{median}}';
    elseif (type==4) % grs, iwue
        wue_data = wue_hourly_summer;
        col_wue = 15;
        wue_name = '{\it WUE_{ei}}';
        units = '(g C hPa kg H_2O^{-1})';
        fig_header = 'iwue_sum_MAH';
        xname = '{\it {MAH(WUE_{ei}}}';
        wue_med_lgd = '{\it WUE_{ei, median}}';
    end
    clear r_agg IGBP_name
    % analyze the top X% windows
    pct_window = [95];
    IGBP_included = unique(IGBP);% % exclude croplands
    idx = find(IGBP_included=='CRO'|IGBP_included=='WET');
    IGBP_included(idx) = [];
    r_Sen_MAH_pct_type = [];
    r_Sen_MAH_pct_agg = [];
    r_Sen_MAH_pct_site = [];
    r_Sen_MAH_pct_base = [];
    trends = [];
    uwue_norm_trend = [];
    mah_sensitivity = [];
    mah_sensitivity_norm = [];
    pft_agg = [];
    site_name_agg = [];
    mah_all_agg = [];
    mah_base_agg = [];
    for i = 1:length(IGBP_included)
        plot_id = [];   
        % for each IGBP
        idx_IGBP = find(IGBP==IGBP_included(i));
        clear site_list
        % declare the matrixs for trend analysis
        sen_agg = nan(length(idx_IGBP), 4, length(pct_window)+2);
        h_agg = nan(length(idx_IGBP), 4, length(pct_window)+2);
        norm_trend_agg = nan(length(idx_IGBP), 4, length(pct_window)+2);
        MAH_trend_agg = nan(length(idx_IGBP), 4);
        for variable = 1:1
            clf;
            annual_change = [];
            count = 0;
            for j = 1:length(idx_IGBP) % loop through sites
                if (site_included(idx_IGBP(j))<166)
                    pct_threshold_type = pct_best_r(i, variable, type);
                    pct_threshold_site = pct_best_r_site(site_included(idx_IGBP(j)), variable, type);
                    pct_threshold_agg = pct_best_r_agg(type);
                    idx = find(wue_data(:,4)==site_included(idx_IGBP(j)));
                    data_tmp = wue_data(idx, :);
                    % define growing season as GPP>0 & Tair>0
                    idx = find(data_tmp(:, 5)>0 & data_tmp(:, 8)>0);
                    grs_data = data_tmp(idx,:);
                    % find the pct values
                    pct_value_type = prctile(grs_data(:,col_wue), pct_threshold_type);
                    pct_value_site = prctile(grs_data(:,col_wue), pct_threshold_site);
                    pct_value_agg = prctile(grs_data(:,col_wue), pct_threshold_agg);
                    pct_value_base = prctile(grs_data(:,col_wue), 50);
                    pct_value_base_low = prctile(grs_data(:,col_wue), 40);
                    pct_value_base_up = prctile(grs_data(:,col_wue), 60);
                    yr_range = unique(grs_data(:, 1));
                    if (length(yr_range)>=10)
                        plot_id = [plot_id; j];
                        count = count+1;
                        clear grs_sum MAH_site MAH_type MAH_agg MAH_base 
                        seasonal_means = [];
                        for yr = 1:length(yr_range) % loop through years
                            idx = find(grs_data(:, 1)==yr_range(yr));
                            data_site_year = grs_data(idx, :);
                            for window = 1:length(pct_window)
                                % use different color to indicate different trends
                                % calculated by different metrics
                                tmp = prctile(data_site_year(:, col_wue),pct_window(window));
                                if (pct_window(window)>=50)
                                    idx = find(data_site_year(:, col_wue)>=tmp);
                                else
                                    idx = find(data_site_year(:, col_wue)<=tmp);
                                end
                                if (variable==2 | variable==4)
                                    grs_sum(yr, window) = sum(data_site_year(idx, col_wue));
                                elseif (variable==1 | variable==3)
                                    grs_sum(yr, window) = median(data_site_year(idx, col_wue));
                                end
                            end
                            if (variable==2 | variable==4)
                                grs_sum(yr, length(pct_window)+1) = sum(data_site_year(:, col_wue),'omitnan');
                            elseif (variable==1 | variable==3)
                                grs_sum(yr, length(pct_window)+1) = median(data_site_year(:, col_wue),'omitnan');
                            end
                            idx = find(data_site_year(:, col_wue)>=pct_value_type);
                            MAH_type(yr) = length(idx);
                            idx = find(data_site_year(:, col_wue)>=pct_value_site);
                            MAH_site(yr) = length(idx);
                            idx = find(data_site_year(:, col_wue)>=pct_value_agg);
                            MAH_agg(yr) = length(idx);
                            idx = find(data_site_year(:, col_wue)>=pct_value_base_low & data_site_year(:, col_wue)<=pct_value_base_up);
                            MAH_base(yr) = length(idx);
                            % output seasonal uwue, MAH_all, MAH_base, gpp, et, and vpd
                            tmp = [median(data_site_year(:, col_wue),'omitnan'), ...
                                length(find(data_site_year(:, col_wue)>=pct_value_agg)),...
                                length(find(data_site_year(:, col_wue)>=pct_value_base_low & data_site_year(:, col_wue)<=pct_value_base_up)),...
                                sum(data_site_year(:, 5),'omitnan'),...
                                sum(data_site_year(:, 6),'omitnan'),...
                                mean(data_site_year(:, 7),'omitnan')];
                            seasonal_means = [seasonal_means; tmp];
                        end
                        annual_change = [annual_change; diff(grs_sum)./mean(grs_sum,'omitnan')];
                    if (variable==1)
                        siteID = dir(site_included(idx_IGBP(j)), 5:10);
                        site_list = {[siteID '(' int2str(length(yr_range)) ')']};
                        IGBP_name{i} = char(IGBP_included(i));
                    end
                        % plot trends and MAH distributions for each site
                        tmp = grs_sum(:, length(pct_window)+1);
                        idx = find(isfinite(tmp)==1);
                        if (length(idx)>=10)
                            tmp = [seasonal_means(idx,:), i*ones(length(idx),1)];
                            mah_sensitivity = [mah_sensitivity; tmp];
                            tmp = [seasonal_means(idx,:)./mean(seasonal_means(idx,:)), ...
                                i*ones(length(idx),1)];
                            mah_sensitivity_norm = [mah_sensitivity_norm; tmp];
                            datain = nan(length(idx),2);
                            datain(:,1) = yr_range(idx);
                            tmp = grs_sum(:, length(pct_window)+1);
                            datain(:,2) = (tmp(idx)-mean(tmp(idx)))/mean(tmp(idx));
                            r = corr(datain(:,2), MAH_type(idx)');
                            r_Sen_MAH_pct_type = [r_Sen_MAH_pct_type; r];
                            r = corr(datain(:,2), MAH_site(idx)');
                            r_Sen_MAH_pct_site = [r_Sen_MAH_pct_site; r];
                            r = corr(datain(:,2), MAH_agg(idx)');
                            r_Sen_MAH_pct_agg = [r_Sen_MAH_pct_agg; r];
                            r = corr(datain(:,2), MAH_base(idx)');
                            r_Sen_MAH_pct_base = [r_Sen_MAH_pct_base; r];
                            site_name_agg = [site_name_agg, site_list];
                            mah_all_agg = [mah_all_agg; MAH_agg(idx)'];
                            mah_base_agg = [mah_base_agg; MAH_base(idx)'];
                            for k = 1:6
                                tmp = [datain(:,1), (seasonal_means(idx, k)-mean(seasonal_means(idx,k)))/mean(seasonal_means(idx,k))];
                                [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] =...
                                    ktaub(tmp, 0.05);
                                trends_tmp(k) = sen*100;
                            end
                            trends = [trends; trends_tmp];
                            [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3 nsigma] =...
                                ktaub(datain, 0.05);
                            uwue_norm_trend = [uwue_norm_trend; sen*100];
                            pft_agg = [pft_agg; i];
                        end
                    end
                end
            end
        end
    end
    % % r(WUE, MAH_all), r(WUE, MAH_base), r(MAH_all, MAH_base)
    clf;
    subplot(2, 2, 1)
    R2 = corr(mah_all_agg, mah_base_agg)^2;
    % one to one
    plot([min(mah_base_agg),max(mah_base_agg)],...
        [min(mah_base_agg),max(mah_base_agg)],'k--')
    hold on
    % line fit
    rg_fit_tmp = polyfit(mah_base_agg, mah_all_agg, 1);
    xxx = linspace(min(mah_base_agg),max(mah_base_agg),10);
    fff = polyval(rg_fit_tmp,xxx);
    plot(xxx,fff,'b-','Linewidth',2)
    ppp = dscatter(mah_base_agg, mah_all_agg);
    colormap('Hot')
    %     c = colorbar('Position', [0.92, 0.17, 0.02, 0.3]);
    axis([min(mah_base_agg) max(mah_base_agg) min(mah_base_agg) max(mah_base_agg)])
    xlabel(['{\it MAH_{uWUE^{40^{th}-60^{th}}}} (h)'], 'FontSize', 12);
    ylabel(['{\it MAH_{uWUE^{78^{th}-99^{th}}}} (h)'], 'FontSize', 12);
    % id = char(erase(id_name(variable+2),'"'));
    t = text(0.02,0.98,['(a) R^2 = ' num2str(sprintf('%1.2f',R2))],...
        'Units', 'Normalized', 'VerticalAlignment', 'Top',...
            'FontSize',12);

    pft = unique(pft_agg);
    idx_plot = [];
    h_type = [];
    tmp = [r_Sen_MAH_pct_base, r_Sen_MAH_pct_agg];
    subplot(2, 2, 2)
    for i = 1:length(pft)
        idx = find(pft_agg==pft(i));
        if (length(idx)>0)
            idx_plot = [idx_plot; i];
        end
        h1 = scatter(r_Sen_MAH_pct_base(idx), r_Sen_MAH_pct_agg(idx), 'o', ...
            'filled', 'MarkerFaceColor', mcolor(i,:));
        hold on
        h_type = [h_type; h1];
        if (i==1)
            xlabel(['Cor({\it MAH_{uWUE^{40^{th}-60^{th}}}}, {\it uWUE_{median}})'], 'FontSize', 12)
            ylabel(['Cor({\it MAH_{uWUE^{78^{th}-99^{th}}}}, {\it uWUE_{median}})'], 'FontSize', 12)
            box on
            plot([min(tmp(:)), max(tmp(:))], [min(tmp(:)), max(tmp(:))], 'k--')
            t = text(0.02,0.98,['(b)'],'Units', 'Normalized',...
                'VerticalAlignment', 'Top','FontSize',12);
        end 
    end
    axis tight
    box on
    subplot(2, 2, 3:4)
    pft = unique(pft_agg);
    idx_plot = [];
    h_type = [];
    tmp = [r_Sen_MAH_pct_agg];
    plot([0 0], [min(tmp) max(tmp)], 'k:')
    hold on
    plot([min(uwue_norm_trend) max(uwue_norm_trend)], [0 0], 'k:')
    axis([min(uwue_norm_trend) max(uwue_norm_trend) min(tmp) max(tmp)])
    for i = 1:length(pft)
        idx = find(pft_agg==pft(i));
        if (length(idx)>0)
            idx_plot = [idx_plot; i];
        end
        h1 = scatter(uwue_norm_trend(idx), tmp(idx), 'o', ...
            'filled', 'MarkerFaceColor', mcolor(i,:));
        hold on
        h_type = [h_type; h1];
        if (i==1)
            h1_tmp = scatter(0, 0, 'o', ...
            'filled', 'MarkerFaceColor', 'k');
            xlabel([wue_med_lgd ' trends ' '(% per year)'], 'FontSize', 12)
            ylabel(['Cor({\it MAH_{uWUE^{78^{th}-99^{th}}}}, {\it uWUE_{median}})'], 'FontSize', 12)
            t = text(0.02,0.98,['(c)'],'Units', 'Normalized',...
                'VerticalAlignment', 'Top','FontSize',12);
        end 
    end
    box on
    ax_copy = axes('Position',get(gca,'Position'),'Visible','Off');
    lgd = legend([h_type],...
                {char(IGBP_name(pft(idx_plot)))},...
                'orientation','horizontal', 'box','off');
    lgd.FontSize = 12;
    lgd.Position = [0.3 0.45 0.3643 0.0452];
    figname=['../plots_hourly/uWUE_trend_full_v7/' ...
                fig_header '_agg_base_comp.jpeg'];
    print('-dpng','-r300',figname);
    % store correlation time series for later use
    if (type==1)
        r_Sen_MAH_pct_agg_uwue = r_Sen_MAH_pct_agg;
        r_Sen_MAH_pct_base_uwue = r_Sen_MAH_pct_base;
    end
end

%% what contribute to the WUE at MAH? how is it different from other percentiles?
% only include the measued data
clc;
close all;
% locate the variables: Tair; Precip; VPD; CO2 mole fraction; GPP; 
% Solar radiation; USTAR
var_test = ["TA_F","P_F","VPD_F",...
    "CO2_F_MDS","SW_IN_F_MDS","USTAR","WS_F","PPFD_IN","TS_F_MDS_1",...
    "SWC_F_MDS_1", "G_F_MDS","NEE_VUT_REF"];
var_idx = [];
for i = 1:length(var_test)
    var_idx = [var_idx; find(strcmpi(var_test(i),var_wanted))];
end
% measured data only
data = [];
for j = 1:length(var_idx)
    col = var_idx(j);
    data(:,j+1) = hourly_data(:,col);
    % QC flags, only use measured data
    idx = find(hourly_data(:,col+1)>0);
    data(idx,j+1) = nan;
end
% set time label
data(:,1) = Time(:,2);
% locate GPP
% umolCO2 m-2 s-1 for hourly data
var_test = ["GPP_NT_VUT_REF"];
var_idx = [];
for i = 1:length(var_test)
    var_idx = [var_idx; find(strcmpi(var_test(i),var_wanted))];
end
% GPP QC
gpp = hourly_data(:,var_idx(1));
% set GPP, umolCO2 m-2 s-1 to gCm-2s-1
data(:,end+1) = gpp*12*(10^-6);
% set ET
data(:,end+1) = daily_ET;
% include site-year info
data(:,end+1) = Time(:,1);
% include site-hour info
data(:,end+1) = Time(:,3);
% include site info
data(:,end+1) = hourly_data(:,end);
vpd_m = [];
Tair_m = [];
cco2_m = [];
% define daytime as solar radiation > 100 W/m2
sw_threshold = 100;
et_threshold = 10^-6;
gpp_threshold = 10^-6;
tair_threshold = 0;
ustar_threshold = 0.1;
uwue_MAH_grs = [];
uwue_MAH_grs_mean = [];
iwue_MAH_grs = [];
iwue_MAH_grs_mean = [];
uwue_MAH_bins = [];
MAH_agg = [];
% declare month DOY
for month = 2:12
    month_start(month) = mon_sum(month-1)+1;
    month_end(month) = mon_sum(month);
end
month_start(1) = 1;
month_end(1) = mon_sum(1);
pct_window = 5:5:95;
clear uwue uwue obs_year pct_value
% calculate WUE based on seasonal means/sums (1 WUE per growing season)
for site = 1:length(site_included)
    lat_tmp = LAT(site);
    % exclude croplands
    if (IGBP(site)~='CRO'|IGBP(site)~='WET')
        count = 1;
        idx = find(data(:,18)==site_included(site));
        data_site = data(idx,:);
        % estimate WUE
        wue = data_site(:,14).*(data_site(:,4).^.5)./data_site(:,15);
        iwue = data_site(:,14).*(data_site(:,4))./data_site(:,15);
        tmp = data_site(:,16);
        tmp(~isfinite(wue)) = [];
        if (length(unique(tmp))>=10)
            yr_range = [min(tmp):max(tmp)];
            yr_start(site) = yr_range(1);
            yr_end(site) = yr_range(end);
            if (length(wue)>0)
                tmp = wue;
                uwue_threshold = [prctile(tmp,pct_best_r_agg(1)),...
                    prctile(tmp,pct_best_r_agg(1)+1),prctile(tmp,50), ...
                    prctile(tmp,40), prctile(tmp,60)];
                for bin = 1:length(pct_window)
                    uwue_threshold_bins(bin) = prctile(tmp,pct_window(bin));
                end
            end
            if (length(iwue)>0)
                tmp = iwue;
                iwue_threshold = [prctile(tmp,pct_best_r_agg(2)),...
                    prctile(tmp,pct_best_r_agg(2)+1),prctile(tmp,50), ...
                    prctile(tmp,40), prctile(tmp,60)];
                for bin = 1:length(pct_window)
                    iwue_threshold_bins(bin) = prctile(tmp,pct_window(bin));
                end
            end
            for yr = 1:length(yr_range)
                idx = find(data_site(:,16)==yr_range(yr));
                data_site_year = data_site(idx,:);
                % precip events
                precip = sum(data_site_year(:,3),'omitnan'); % total precip [mm]
                % precip events
                idx = find(data_site_year(:,3)>0);
                doy_tmp = unique(data_site_year(idx,1));
                % filter out data during and one day after a precip event
                for dd = 1:length(doy_tmp)
                    idx = find(data_site_year(:,1)==doy_tmp(dd)|...
                        data_site_year(:,1)==(doy_tmp(dd)+1));
                    data_site_year(idx,:) = nan;
                end
                % filter out low ustar
                idx = find(data_site_year(:,7)<ustar_threshold);
                data_site_year(idx,:) = nan;
                % filter out nighttime data
                idx = find(data_site_year(:,6)<sw_threshold);
                data_site_year(idx,:) = nan;
                % condensation events
                idx = find(data_site_year(:,15)<et_threshold);
                data_site_year(idx,:) = nan;
                % set growing season as the period when daily GPP>5%
                % daily GPP max and daily Tair>0
                idx = isfinite(data_site_year(:, 1));
                days = unique(data_site_year(idx, 1));
                tmp = [];
                for day = 1:length(days)
                    idx = find(data_site_year(:, 1)==days(day));
                    if (length(idx)>0)
                        % DOY, Tair, GPP
                        xxx = [data_site_year(idx(1), 1), ...
                            mean(data_site_year(idx,2), 'omitnan'),...
                            sum(data_site_year(idx,14), 'omitnan')];
                        tmp = [tmp; xxx];
                    end
                end
                idx = find(tmp(:,2)<=0 | tmp(:,3)<max(tmp(:,3))*0.05);
                if (length(idx)>1)
                    days = unique(tmp(idx, 1));
                    for day = 1:length(days)
                        idx = find(data_site_year(:, 1)==days(day));
                        if (length(idx)>0)
                            data_site_year(idx,:) = nan;
                        end
                    end
                end
                % When GPP<0
                idx = find(data_site_year(:,14)<gpp_threshold);
                data_site_year(idx,:) = nan;
                % When Tair<0
                idx = find(data_site_year(:,2)<tair_threshold);
                data_site_year(idx,:) = nan;
                uwue_hourly = data_site_year(:,14).*(data_site_year(:,4).^.5)./data_site_year(:,15);
                iwue_hourly = data_site_year(:,14).*(data_site_year(:,4))./data_site_year(:,15);
                % daily aggregate
                % output year, DOY, WUE, site, GPP, ET, VPD, Tair, [CO2], USTAR
                idx = isfinite(data_site_year(:, 1));
                days = unique(data_site_year(idx, 1));
                tmp = [];
                for day = 1:length(days)
                    idx = find(data_site_year(:, 1)==days(day));
                    if (length(idx)>0)
                        xxx = [data_site_year(idx(1),10), data_site_year(idx(1), 1), ...
                        mean(uwue_hourly(idx), 'omitnan'), data_site_year(idx(1),12), ...
                        sum(data_site_year(idx,14), 'omitnan'), sum(data_site_year(idx,15), 'omitnan'), ...
                        mean(data_site_year(idx,4), 'omitnan'), mean(data_site_year(idx,2), 'omitnan'), ...
                        mean(data_site_year(idx,5), 'omitnan'), mean(data_site_year(idx,7), 'omitnan')];
                        xxx(3) = xxx(5)*(xxx(7)^.5)/xxx(6); % uWUE from daily agg values
                        xxx(end+1) = xxx(5)*(xxx(7))/xxx(6); % WUEei from daily agg values
                        tmp = [tmp; xxx];
                    end
                end
                uwue_daily = tmp;
                if (length(uwue_daily)>0) % data exist
                    % output year, DOY, WUE, site, GPP, ET, VPD, Tair, [CO2], USTAR
                    if (sum(isfinite(uwue_daily(:,3)))>1) 
                        % growing-season DOY
                        doy = data_site_year(:,1);
                        doy(~isfinite(doy)) = [];
%                         % uWUE
                        tmp = [uwue_hourly, data_site_year, data_site_year(:,[2:15])-mean(data_site_year(:,[2:15]),'omitnan')];
                        % growing season mean with number of MAH
                        idx = find(uwue_hourly>=uwue_threshold(1));
                        MAH = length(idx);
                        % MAHbase (uwue40<uwue<uwue60)
                        idx = find(uwue_hourly>=uwue_threshold(4)&uwue_hourly<=uwue_threshold(5));
                        baseline_hour = length(idx);                        
                        % MAHtran (>50th WUE & < pct_best_r_agg)
                        idx = find(uwue_hourly>=uwue_threshold(3)&uwue_hourly<uwue_threshold(1));
                        tran_hour = length(idx);
                        % GRS WUE from aggregated GPP, ET and VPD
                        grs_wue = sum(data_site_year(:,14), 'omitnan')*(mean(data_site_year(:,4), 'omitnan')^.5)...
                            /sum(data_site_year(:,15), 'omitnan');
                        output = [mean(tmp, 'omitnan'), grs_wue, MAH, ...
                            baseline_hour, tran_hour, length(unique(doy))];
                        % modify the columns in outputs
                        % median uWUE, precip sum, GPP sum, ET sum
                        output(:,1) = median(uwue_hourly, 'omitnan');
                        output(:,4) = precip;
                        output(:,15) = sum(data_site_year(:,14), 'omitnan');
                        output(:,16) = sum(data_site_year(:,15), 'omitnan');
                        uwue_MAH_grs_mean = [uwue_MAH_grs_mean; output];
                        % WUEei
                        tmp = [iwue_hourly, data_site_year, data_site_year(:,[2:15])-mean(data_site_year(:,[2:15]),'omitnan')];
                        % growing season mean with number of MAH
                        idx = find(iwue_hourly>=iwue_threshold(1));
                        MAH = length(idx);
                        % MAHbase (uwue40<uwue<uwue60)
                        idx = find(iwue_hourly>=iwue_threshold(4)&iwue_hourly<=iwue_threshold(5));
                        baseline_hour = length(idx);                        
                        % MAHtran (>50th WUE & < pct_best_r_agg)
                        idx = find(iwue_hourly>=iwue_threshold(3)&iwue_hourly<iwue_threshold(1));
                        tran_hour = length(idx);
                        % GRS WUE from aggregated GPP, ET and VPD
                        grs_wue = sum(data_site_year(:,14), 'omitnan')*(mean(data_site_year(:,4), 'omitnan'))...
                            /sum(data_site_year(:,15), 'omitnan');
                        output = [mean(tmp, 'omitnan'), grs_wue, MAH, ...
                            baseline_hour, tran_hour, length(unique(doy))];
                        % modify the columns in outputs
                        % median uWUE, precip sum, GPP sum, ET sum
                        output(:,1) = median(iwue_hourly, 'omitnan');
                        output(:,4) = precip;
                        output(:,15) = sum(data_site_year(:,14), 'omitnan');
                        output(:,16) = sum(data_site_year(:,15), 'omitnan');
                        iwue_MAH_grs_mean = [iwue_MAH_grs_mean; output];
                    end
                end
            end
        end
    end
end

%% what contribute to the WUE at MAH? how is it different from other percentiles?
% summer time only
% only include the measued data
clc;
close all;
% locate the variables: Tair; Precip; VPD; CO2 mole fraction; GPP; 
% Solar radiation; USTAR
var_test = ["TA_F","P_F","VPD_F",...
    "CO2_F_MDS","SW_IN_F_MDS","USTAR","WS_F","PPFD_IN","TS_F_MDS_1",...
    "SWC_F_MDS_1", "G_F_MDS","NEE_VUT_REF"];
var_idx = [];
for i = 1:length(var_test)
    var_idx = [var_idx; find(strcmpi(var_test(i),var_wanted))];
end
% measured data only
data = [];
for j = 1:length(var_idx)
    col = var_idx(j);
    data(:,j+1) = hourly_data(:,col);
    % QC flags, only use measured data
    idx = find(hourly_data(:,col+1)>0);
    data(idx,j+1) = nan;
end
% set time label
data(:,1) = Time(:,2);
% locate GPP
% umolCO2 m-2 s-1 for hourly data
var_test = ["GPP_NT_VUT_REF"];
var_idx = [];
for i = 1:length(var_test)
    var_idx = [var_idx; find(strcmpi(var_test(i),var_wanted))];
end
% GPP QC
gpp = hourly_data(:,var_idx(1));
% set GPP, umolCO2 m-2 s-1 to gCm-2s-1
data(:,end+1) = gpp*12*(10^-6);
% set ET
data(:,end+1) = daily_ET;
% include site-year info
data(:,end+1) = Time(:,1);
% include site-hour info
data(:,end+1) = Time(:,3);
% include site info
data(:,end+1) = hourly_data(:,end);
vpd_m = [];
Tair_m = [];
cco2_m = [];
% define daytime as solar radiation > 100 W/m2
sw_threshold = 100;
et_threshold = 10^-6;
gpp_threshold = 10^-6;
tair_threshold = 0;
ustar_threshold = 0.1;
uwue_MAH_sum_mean = [];
iwue_MAH_sum_mean = [];
% declare month DOY
for month = 2:12
    month_start(month) = mon_sum(month-1)+1;
    month_end(month) = mon_sum(month);
end
month_start(1) = 1;
month_end(1) = mon_sum(1);
pct_window = 5:5:95;
clear uwue uwue obs_year pct_value
% calculate WUE based on seasonal means/sums (1 WUE per growing season)
for site = 1:length(site_included)
    lat_tmp = LAT(site);
    % exclude croplands
    if (IGBP(site)~='CRO'|IGBP(site)~='WET')
        count = 1;
        idx = find(data(:,18)==site_included(site));
        data_site = data(idx,:);
        % estimate WUE
        wue = data_site(:,14).*(data_site(:,4).^.5)./data_site(:,15);
        iwue = data_site(:,14).*(data_site(:,4))./data_site(:,15);
        tmp = data_site(:,16);
        tmp(~isfinite(wue)) = [];
        if (length(unique(tmp))>=10)
            yr_range = [min(tmp):max(tmp)];
            yr_start(site) = yr_range(1);
            yr_end(site) = yr_range(end);
            % WUE for the entire time series
            if (length(wue)>0)
                tmp = wue;
                uwue_threshold = [prctile(tmp,pct_best_r_agg(1)),...
                    prctile(tmp,pct_best_r_agg(1)+1),prctile(tmp,50), ...
                    prctile(tmp,40), prctile(tmp,60)];
                for bin = 1:length(pct_window)
                    uwue_threshold_bins(bin) = prctile(tmp,pct_window(bin));
                end
            end
            if (length(iwue)>0)
                tmp = iwue;
                iwue_threshold = [prctile(tmp,pct_best_r_agg(2)),...
                    prctile(tmp,pct_best_r_agg(2)+1),prctile(tmp,50), ...
                    prctile(tmp,40), prctile(tmp,60)];
                for bin = 1:length(pct_window)
                    iwue_threshold_bins(bin) = prctile(tmp,pct_window(bin));
                end
            end
            for yr = 1:length(yr_range)
                idx = find(data_site(:,16)==yr_range(yr));
                data_site_year = data_site(idx,:);
                % precip events
                precip = sum(data_site_year(:,3),'omitnan'); % total precip [mm]
                % precip events
                idx = find(data_site_year(:,3)>0);
                doy_tmp = unique(data_site_year(idx,1));
                % filter out data during and one day after a precip event
                for dd = 1:length(doy_tmp)
                    idx = find(data_site_year(:,1)==doy_tmp(dd)|...
                        data_site_year(:,1)==(doy_tmp(dd)+1));
                    data_site_year(idx,:) = nan;
                end
                % filter out low ustar
                idx = find(data_site_year(:,7)<ustar_threshold);
                data_site_year(idx,:) = nan;
                % filter out nighttime data
                idx = find(data_site_year(:,6)<sw_threshold);
                data_site_year(idx,:) = nan;
                % condensation events
                idx = find(data_site_year(:,15)<et_threshold);
                data_site_year(idx,:) = nan;
                % define summer season
                if (lat_tmp>0)
                    summer_start = mon_sum(5)+1;
                    summer_end = mon_sum(8);
                    idx_other = find(data_site_year(:, 1)<summer_start | ...
                        data_site_year(:, 1)>summer_end);
                else
                    summer_start = mon_sum(11)+1;
                    summer_end = mon_sum(2);
                    idx_other = find(data_site_year(:, 1)<summer_start & ...
                        data_site_year(:, 1)>summer_end);
                end
                data_site_year(idx_other,:) = nan;
                % set growing season as the period when daily GPP>5%
                % daily GPP max and daily Tair>0
                idx = isfinite(data_site_year(:, 1));
                days = unique(data_site_year(idx, 1));
                tmp = [];
                for day = 1:length(days)
                    idx = find(data_site_year(:, 1)==days(day));
                    if (length(idx)>0)
                        % DOY, Tair, GPP
                        xxx = [data_site_year(idx(1), 1), ...
                            mean(data_site_year(idx,2), 'omitnan'),...
                            sum(data_site_year(idx,14), 'omitnan')];
                        tmp = [tmp; xxx];
                    end
                end
                idx = find(tmp(:,2)<=0 | tmp(:,3)<max(tmp(:,3))*0.05);
                if (length(idx)>1)
                    days = unique(tmp(idx, 1));
                    for day = 1:length(days)
                        idx = find(data_site_year(:, 1)==days(day));
                        if (length(idx)>0)
                            data_site_year(idx,:) = nan;
                        end
                    end
                end
                % When GPP<0
                idx = find(data_site_year(:,14)<gpp_threshold);
                data_site_year(idx,:) = nan;
                % When Tair<0
                idx = find(data_site_year(:,2)<tair_threshold);
                data_site_year(idx,:) = nan;
                uwue_hourly = data_site_year(:,14).*(data_site_year(:,4).^.5)./data_site_year(:,15);
                iwue_hourly = data_site_year(:,14).*(data_site_year(:,4))./data_site_year(:,15);
                % daily aggregate
                % output year, DOY, WUE, site, GPP, ET, VPD, Tair, [CO2], USTAR
                idx = isfinite(data_site_year(:, 1));
                days = unique(data_site_year(idx, 1));
                tmp = [];
                for day = 1:length(days)
                    idx = find(data_site_year(:, 1)==days(day));
                    if (length(idx)>0)
                        xxx = [data_site_year(idx(1),10), data_site_year(idx(1), 1), ...
                        mean(uwue_hourly(idx), 'omitnan'), data_site_year(idx(1),12), ...
                        sum(data_site_year(idx,14), 'omitnan'), sum(data_site_year(idx,15), 'omitnan'), ...
                        mean(data_site_year(idx,4), 'omitnan'), mean(data_site_year(idx,2), 'omitnan'), ...
                        mean(data_site_year(idx,5), 'omitnan'), mean(data_site_year(idx,7), 'omitnan')];
                        xxx(3) = xxx(5)*(xxx(7)^.5)/xxx(6); % uWUE from daily agg values
                        xxx(end+1) = xxx(5)*(xxx(7))/xxx(6); % WUEei from daily agg values
                        tmp = [tmp; xxx];
                    end
                end
                uwue_daily = tmp;
                if (length(uwue_daily)>0) % data exist
                    % output year, DOY, WUE, site, GPP, ET, VPD, Tair, [CO2], USTAR
                    if (sum(isfinite(uwue_daily(:,3)))>1) 
                        % growing-season DOY
                        doy = data_site_year(:,1);
                        doy(~isfinite(doy)) = [];
                        % uWUE
                        tmp = [uwue_hourly, data_site_year, data_site_year(:,[2:15])-mean(data_site_year(:,[2:15]),'omitnan')];
                        % growing season mean with number of MAH
                        idx = find(uwue_hourly>=uwue_threshold(1));
                        MAH = length(idx);
                        % MAHbase (uwue40<uwue<uwue60)
                        idx = find(uwue_hourly>=uwue_threshold(4)&uwue_hourly<=uwue_threshold(5));
                        baseline_hour = length(idx);                        
                        % MAHtran (>50th WUE & < pct_best_r_agg)
                        idx = find(uwue_hourly>=uwue_threshold(3)&uwue_hourly<uwue_threshold(1));
                        tran_hour = length(idx);
                        % GRS WUE from aggregated GPP, ET and VPD
                        grs_wue = sum(data_site_year(:,14), 'omitnan')*(mean(data_site_year(:,4), 'omitnan')^.5)...
                            /sum(data_site_year(:,15), 'omitnan');
                        output = [mean(tmp, 'omitnan'), grs_wue, MAH, ...
                            baseline_hour, tran_hour, length(unique(doy))];
                        % modify the columns in outputs
                        % median uWUE, precip sum, GPP sum, ET sum
                        output(:,1) = median(uwue_hourly, 'omitnan');
                        output(:,4) = precip;
                        output(:,15) = sum(data_site_year(:,14), 'omitnan');
                        output(:,16) = sum(data_site_year(:,15), 'omitnan');
                        uwue_MAH_sum_mean = [uwue_MAH_sum_mean; output];
                        % WUEei
                        tmp = [iwue_hourly, data_site_year, data_site_year(:,[2:15])-mean(data_site_year(:,[2:15]),'omitnan')];
                        % growing season mean with number of MAH
                        idx = find(iwue_hourly>=iwue_threshold(1));
                        MAH = length(idx);
                        % MAHbase (uwue40<uwue<uwue60)
                        idx = find(iwue_hourly>=iwue_threshold(4)&iwue_hourly<=iwue_threshold(5));
                        baseline_hour = length(idx);                        
                        % MAHtran (>50th WUE & < pct_best_r_agg)
                        idx = find(iwue_hourly>=iwue_threshold(3)&iwue_hourly<iwue_threshold(1));
                        tran_hour = length(idx);
                        % GRS WUE from aggregated GPP, ET and VPD
                        grs_wue = sum(data_site_year(:,14), 'omitnan')*(mean(data_site_year(:,4), 'omitnan'))...
                            /sum(data_site_year(:,15), 'omitnan');
                        output = [mean(tmp, 'omitnan'), grs_wue, MAH, ...
                            baseline_hour, tran_hour, length(unique(doy))];
                        % modify the columns in outputs
                        % median uWUE, precip sum, GPP sum, ET sum
                        output(:,1) = median(iwue_hourly, 'omitnan');
                        output(:,4) = precip;
                        output(:,15) = sum(data_site_year(:,14), 'omitnan');
                        output(:,16) = sum(data_site_year(:,15), 'omitnan');
                        iwue_MAH_sum_mean = [iwue_MAH_sum_mean; output];
                    end
                end
            end
        end
    end
end

%% analyze results when there are moren than 10 years of measurements
clc;
IGBP_included = unique(IGBP);% % exclude croplands
idx = find(IGBP_included=='CRO'|IGBP_included=='WET');
IGBP_included(idx) = [];
uwue_MAH_grs_mean_subset = [];
iwue_MAH_grs_mean_subset = [];
uwue_MAH_sum_mean_subset = [];
iwue_MAH_sum_mean_subset = [];
for i = 1:length(IGBP_included)
    % find sites within an ecosystem type
    idx_IGBP = find(IGBP==IGBP_included(i));
    for j = 1:length(idx_IGBP)
        site = site_included(idx_IGBP(j));
        % uwue_MAH_grs_mean
        idx = find(uwue_MAH_grs_mean(:,19)==site);
        site_year = unique(uwue_MAH_grs_mean(idx,17));
        % include the site data if there's more than 10 site-years
        if (length(site_year)>=10)
            tmp = [uwue_MAH_grs_mean(idx,:), ones(length(idx),1)*i]; % attach type info
            uwue_MAH_grs_mean_subset = [uwue_MAH_grs_mean_subset; tmp];
        end
        % iwue_MAH_grs_mean
        idx = find(iwue_MAH_grs_mean(:,19)==site);
        site_year = unique(iwue_MAH_grs_mean(idx,17));
        % include the site data if there's more than 10 site-years
        if (length(site_year)>=10)
            tmp = [iwue_MAH_grs_mean(idx,:), ones(length(idx),1)*i]; % attach type info
            iwue_MAH_grs_mean_subset = [iwue_MAH_grs_mean_subset; tmp];
        end
        % uwue_MAH_sum_mean
        idx = find(uwue_MAH_sum_mean(:,19)==site);
        site_year = unique(uwue_MAH_sum_mean(idx,17));
        % include the site data if there's more than 10 site-years
        if (length(site_year)>=10)
            tmp = [uwue_MAH_sum_mean(idx,:), ones(length(idx),1)*i]; % attach type info
            uwue_MAH_sum_mean_subset = [uwue_MAH_sum_mean_subset; tmp];
        end
        % iwue_MAH_sum_mean
        idx = find(iwue_MAH_sum_mean(:,19)==site);
        site_year = unique(iwue_MAH_sum_mean(idx,17));
        % include the site data if there's more than 10 site-years
        if (length(site_year)>=10)
            tmp = [iwue_MAH_sum_mean(idx,:), ones(length(idx),1)*i]; % attach type info
            iwue_MAH_sum_mean_subset = [iwue_MAH_sum_mean_subset; tmp];
        end
    end
end
%% random-forest model variable importance: which factor controls grs-mean WUE?
clc;
close all
units = '(g C hPa^{0.5} kg H_2O^{-1})';
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", ...
    "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", ...
    "(q)", "(r)", "(s)", "(t)", "(u)", "(v)", "(w)", "(x)"];
for type = 1:4
    clf;
    if (type==1)
        wue_data = uwue_MAH_grs_mean_subset;
        wue_name = '{\it uWUE}';
        wue_med_lgd = '{\it uWUE_{median}}';
        wue_agg_lgd = '{\it uWUE_{aggregate}}';
        fig_header = 'RF_VR_WUE_ctrl_uWUE_grs';
    elseif (type==2)
        wue_data = iwue_MAH_grs_mean_subset;
        wue_name = '{\it WUE_{ei}}';
        wue_med_lgd = '{\it WUE_{ei, median}}';
        wue_agg_lgd = '{\it WUE_{ei, aggregate}}';
        fig_header = 'RF_VR_WUE_ctrl_iWUE_grs';
    elseif (type==3)
        wue_data = uwue_MAH_sum_mean_subset;
        wue_name = '{\it uWUE}';
        wue_med_lgd = '{\it uWUE_{median}}';
        wue_agg_lgd = '{\it uWUE_{aggregate}}';
        fig_header = 'RF_VR_WUE_ctrl_uWUE_sum';
    elseif (type==4)
        wue_data = iwue_MAH_sum_mean_subset;
        wue_name = '{\it WUE_{ei}}';
        wue_med_lgd = '{\it WUE_{ei, median}}';
        wue_agg_lgd = '{\it WUE_{ei, aggregate}}';
        fig_header = 'RF_VR_WUE_ctrl_iWUE_sum';
    end
    for variable = 1:2
        if (variable==1)
            yname = wue_med_lgd;
            WUE = wue_data(:,1);
        elseif (variable==2)
            yname = wue_agg_lgd;
            WUE = wue_data(:,34);
        end
        idx = find(isfinite(WUE)==1);
        WUE = WUE(idx);
        Tair = wue_data(idx,3);
        VPD = wue_data(idx,5);
        CO2 = wue_data(idx,6);
        SW = wue_data(idx,7);
        USTAR = wue_data(idx,8);
        Wind = wue_data(idx,9);
        PPFD = wue_data(idx,10);
        Tsoil = wue_data(idx,11);
        SWC = wue_data(idx,12);
        G = wue_data(idx,13);
        NEE = wue_data(idx,14);
        GPP = wue_data(idx,15);
        ET = wue_data(idx,16);
        MAH_all = wue_data(idx,35);
        MAH_base = wue_data(idx,36);
        obs_length = wue_data(idx,38);
        DOY = categorical(round(wue_data(idx,2)));
        Hour = categorical(round(wue_data(idx,18)));
        Site = categorical(wue_data(idx,19));
        Type = categorical(wue_data(idx,39));
        X = table(Tair, VPD, CO2, SW, USTAR, Tsoil, SWC, G,...
            GPP, ET, Site, Type, ....
            MAH_all, MAH_base, obs_length,...
            WUE);
        t = templateTree('NumPredictorsToSample','all', ...
            'PredictorSelection', 'curvature', 'Surrogate', 'on');
        rng('default'); % For reproducibility
        Mdl = fitrensemble(X, 'WUE', 'Method', 'bag', 'NumLearningCycles', 100, 'Learners', t);
        imp = oobPermutedPredictorImportance(Mdl);
        subplot(2,2,variable)
        [tmp, idx] = sort(imp,'descend');
        col = length(tmp);
        bar(tmp(1:col));
        xlabel('Predictors', 'FontSize', 12);
        % adjust xtick format
        var_name = Mdl.PredictorNames;
        var_name{5} = '{\it u*}';
        var_name{13} = '{\it MAH_{{78-99}}}';
        var_name{14} = '{\it MAH_{{40-60}}}';
        var_name{15} = '{Obs length}';
        set(gca,'xtick',[1:col])
        set(gca,'xticklabel',string(var_name(idx)), 'XTickLabelRotation', -45)
        ylim([0 max(imp)+2])
        ylabel('Predictor importance estimates', 'FontSize', 12);
        id = char(erase(id_name(variable),'"'));
        t = text(0.02,0.98,[id ' ' yname],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12);
        subplot(2,2,variable+2)
        % R2
        yHat = oobPredict(Mdl);
        R2 = corr(Mdl.Y,yHat)^2;
        % one to one
        plot([min(Mdl.Y),max(Mdl.Y)],[min(Mdl.Y),max(Mdl.Y)],'k--')
        hold on
        % line fit
        rg_fit_tmp = polyfit(Mdl.Y, yHat, 1);
        xxx = linspace(min(Mdl.Y),max(Mdl.Y),10);
        fff = polyval(rg_fit_tmp,xxx);
        plot(xxx,fff,'b-','Linewidth',2)
        ppp = dscatter(Mdl.Y, yHat);
        colormap('Hot')
        axis([min(Mdl.Y) max(Mdl.Y) min(Mdl.Y) max(Mdl.Y)])
        xlabel(['Observed ' yname], 'FontSize', 12);
        ylabel(['Estimated ' yname], 'FontSize', 12);
        id = char(erase(id_name(variable+2),'"'));
        t = text(0.02,0.98,[id ' ' yname],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                'FontSize',12);
        t = text(0.02,0.78,['R^2 = ' num2str(sprintf('%1.2f',R2))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                'FontSize',12);
    end
    figname=['../plots_hourly/uWUE_trend_full_v7/' fig_header '.jpeg'];
    print('-djpeg','-r300',[figname]);
end

%% random-forest model variable importance: which factor controls grs-mean WUE IAV?
clc;
close all
units = '(g C hPa^{0.5} kg H_2O^{-1})';
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", ...
    "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", ...
    "(q)", "(r)", "(s)", "(t)", "(u)", "(v)", "(w)", "(x)"];
for type = 1:4
    clf;
    if (type==1)
        wue_data = uwue_MAH_grs_mean_subset;
        wue_name = '{\it uWUE}';
        wue_med_lgd = '{\it uWUE_{median}}';
        wue_agg_lgd = '{\it uWUE_{aggregate}}';
        fig_header = 'RF_VR_WUE_IAV_uWUE_grs';
    elseif (type==2)
        wue_data = iwue_MAH_grs_mean_subset;
        wue_name = '{\it WUE_{ei}}';
        wue_med_lgd = '{\it WUE_{ei, median}}';
        wue_agg_lgd = '{\it WUE_{ei, aggregate}}';
        fig_header = 'RF_VR_WUE_IAV_iWUE_grs';
    elseif (type==3)
        wue_data = uwue_MAH_sum_mean_subset;
        wue_name = '{\it uWUE}';
        wue_med_lgd = '{\it uWUE_{median}}';
        wue_agg_lgd = '{\it uWUE_{aggregate}}';
        fig_header = 'RF_VR_WUE_IAV_uWUE_sum';
    elseif (type==4)
        wue_data = iwue_MAH_sum_mean_subset;
        wue_name = '{\it WUE_{ei}}';
        wue_med_lgd = '{\it WUE_{ei, median}}';
        wue_agg_lgd = '{\it WUE_{ei, aggregate}}';
        fig_header = 'RF_VR_WUE_IAV_iWUE_sum';
    end
    % normalized trend estimates
    grs_trend_agg = [];
    for i = 1:length(site_included)
        idx = find(wue_data(:, 19)==site_included(i));
        data_site = wue_data(idx,:);
        tmp = (data_site(2:end,:)-data_site(1:end-1,:))./...
            (data_site(2:end,17)-data_site(1:end-1,17))./...
            mean(data_site,'omitnan');
        % include site and type info
        if (length(tmp)>0)
            tmp(:,19) = site_included(i);
            tmp(:,39) = data_site(1, 39);
        end
        grs_trend_agg = [grs_trend_agg; tmp];
    end
    % convert to % per year
    grs_trend_agg = grs_trend_agg*100;
    for variable = 1:2
        if (variable==1)
            yname = wue_med_lgd;
            WUE = grs_trend_agg(:,1);
        elseif (variable==2)
            yname = wue_agg_lgd;
            WUE = grs_trend_agg(:,34);
        end
        idx = find(isfinite(WUE)==1);
        WUE = WUE(idx);
        Tair = grs_trend_agg(idx,3);
        VPD = grs_trend_agg(idx,5);
        CO2 = grs_trend_agg(idx,6);
        SW = grs_trend_agg(idx,7);
        USTAR = grs_trend_agg(idx,8);
        Wind = grs_trend_agg(idx,9);
        PPFD = grs_trend_agg(idx,10);
        Tsoil = grs_trend_agg(idx,11);
        SWC = grs_trend_agg(idx,12);
        G = grs_trend_agg(idx,13);
        NEE = grs_trend_agg(idx,14);
        GPP = grs_trend_agg(idx,15);
        ET = grs_trend_agg(idx,16);
        MAH_all = grs_trend_agg(idx,35);
        MAH_base = grs_trend_agg(idx,36);
        obs_length = grs_trend_agg(idx,38);
        DOY = categorical(round(grs_trend_agg(idx,2)));
        Hour = categorical(round(grs_trend_agg(idx,18)));
        Site = categorical(grs_trend_agg(idx,19));
        Type = categorical(grs_trend_agg(idx,39));
        X = table(Tair, VPD, CO2, SW, USTAR, Tsoil, SWC, G,...
            GPP, ET, Site, Type, ....
            MAH_all, MAH_base, obs_length,...
            WUE);
        t = templateTree('NumPredictorsToSample','all', ...
            'PredictorSelection', 'curvature', 'Surrogate', 'on');
        rng('default'); % For reproducibility
        Mdl = fitrensemble(X, 'WUE', 'Method', 'bag', 'NumLearningCycles', 100, 'Learners', t);
        imp = oobPermutedPredictorImportance(Mdl);
        subplot(2,2,variable)
        [tmp, idx] = sort(imp,'descend');
        col = length(tmp);
        bar(tmp(1:col));
        xlabel('Predictors', 'FontSize', 12);
        % adjust xtick format
        var_name = Mdl.PredictorNames;
        var_name{5} = '{\it u*}';
        var_name{13} = '{\it MAH_{{78-99}}}';
        var_name{14} = '{\it MAH_{{40-60}}}';
        var_name{15} = '{Obs length}';
        set(gca,'xtick',[1:col])
        set(gca,'xticklabel',string(var_name(idx)), 'XTickLabelRotation', -45)
        ylim([0 max(imp)+1])
        ylabel('Predictor importance estimates', 'FontSize', 12);
        xlabel('Normalized IAV (% yr^{-1})', 'FontSize', 12);
        id = char(erase(id_name(variable),'"'));
        t = text(0.02,0.98,[id ' ' yname ' IAV'],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12);
        subplot(2,2,variable+2)
        % R2
        yHat = oobPredict(Mdl);
        R2 = corr(Mdl.Y,yHat)^2;
        % one to one
        plot([min(Mdl.Y),max(Mdl.Y)],[min(Mdl.Y),max(Mdl.Y)],'k--')
        hold on
        % line fit
        rg_fit_tmp = polyfit(Mdl.Y, yHat, 1);
        xxx = linspace(min(Mdl.Y),max(Mdl.Y),10);
        fff = polyval(rg_fit_tmp,xxx);
        plot(xxx,fff,'b-','Linewidth',2)
        ppp = dscatter(Mdl.Y, yHat);
        colormap('Hot')
        axis([min(Mdl.Y) max(Mdl.Y) min(Mdl.Y) max(Mdl.Y)])
        xlabel(['Observed ' yname ' IAV (% yr^{-1})'], 'FontSize', 12);
        ylabel(['Estimated ' yname ' IAV (% yr^{-1})'], 'FontSize', 12);
        id = char(erase(id_name(variable+2),'"'));
        t = text(0.02,0.98,[id ' ' yname],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                'FontSize',12);
        t = text(0.02,0.78,['R^2 = ' num2str(sprintf('%1.2f',R2))],...
            'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                'FontSize',12);
    end
    figname=['../plots_hourly/uWUE_trend_full_v7/' fig_header '.jpeg'];
    print('-djpeg','-r300',[figname]);
end

%% random-forest model variable importance: which factor controls grs-mean WUE IAV?
clc;
close all
units = '(g C hPa^{0.5} kg H_2O^{-1})';
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", ...
    "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", ...
    "(q)", "(r)", "(s)", "(t)", "(u)", "(v)", "(w)", "(x)"];
for type = 1:4
    clf;
    if (type==1)
        wue_data = uwue_MAH_grs_mean_subset;
        wue_name = '{\it uWUE}';
        wue_med_lgd = '{\it uWUE_{median}}';
        wue_agg_lgd = '{\it uWUE_{aggregate}}';
        fig_header = 'RF_VR_WUE_IAV_uWUE_grs';
    elseif (type==2)
        wue_data = iwue_MAH_grs_mean_subset;
        wue_name = '{\it WUE_{ei}}';
        wue_med_lgd = '{\it WUE_{ei, median}}';
        wue_agg_lgd = '{\it WUE_{ei, aggregate}}';
        fig_header = 'RF_VR_WUE_IAV_iWUE_grs';
    elseif (type==3)
        wue_data = uwue_MAH_sum_mean_subset;
        wue_name = '{\it uWUE}';
        wue_med_lgd = '{\it uWUE_{median}}';
        wue_agg_lgd = '{\it uWUE_{aggregate}}';
        fig_header = 'RF_VR_WUE_IAV_uWUE_sum';
    elseif (type==4)
        wue_data = iwue_MAH_sum_mean_subset;
        wue_name = '{\it WUE_{ei}}';
        wue_med_lgd = '{\it WUE_{ei, median}}';
        wue_agg_lgd = '{\it WUE_{ei, aggregate}}';
        fig_header = 'RF_VR_WUE_IAV_iWUE_sum';
    end
    % normalized trend estimates
    grs_trend_agg = [];
    for i = 1:length(site_included)
        idx = find(wue_data(:, 19)==site_included(i));
        data_site = wue_data(idx,:);
        tmp = (data_site(2:end,:)-data_site(1:end-1,:))./...
            (data_site(2:end,17)-data_site(1:end-1,17))./...
            mean(data_site,'omitnan');
        % include site and type info
        if (length(tmp)>0)
            tmp(:,19) = site_included(i);
            tmp(:,39) = data_site(1, 39);
        end
        grs_trend_agg = [grs_trend_agg; tmp];
    end
    % convert to % per year
    grs_trend_agg = grs_trend_agg*100;
    for mah = 1:3
        clf;
        for variable = 1:2
            if (variable==1)
                yname = wue_med_lgd;
                WUE = grs_trend_agg(:,1);
            elseif (variable==2)
                yname = wue_agg_lgd;
                WUE = grs_trend_agg(:,34);
            end
            idx = find(isfinite(WUE)==1);
            WUE = WUE(idx);
            Tair = grs_trend_agg(idx,3);
            VPD = grs_trend_agg(idx,5);
            CO2 = grs_trend_agg(idx,6);
            SW = grs_trend_agg(idx,7);
            USTAR = grs_trend_agg(idx,8);
            Wind = grs_trend_agg(idx,9);
            PPFD = grs_trend_agg(idx,10);
            Tsoil = grs_trend_agg(idx,11);
            SWC = grs_trend_agg(idx,12);
            G = grs_trend_agg(idx,13);
            NEE = grs_trend_agg(idx,14);
            GPP = grs_trend_agg(idx,15);
            ET = grs_trend_agg(idx,16);
            MAH_all = grs_trend_agg(idx,35);
            MAH_base = grs_trend_agg(idx,36);
            MAH_tran = grs_trend_agg(idx,37);
            obs_length = grs_trend_agg(idx,38);
            DOY = categorical(round(grs_trend_agg(idx,2)));
            Hour = categorical(round(grs_trend_agg(idx,18)));
            Site = categorical(grs_trend_agg(idx,19));
            Type = categorical(grs_trend_agg(idx,39));
            % evaluate different MAH difinition
            if (mah==1)
                fig_desp=['_MAHall'];
                X = table(Tair, VPD, CO2, SW, USTAR, Tsoil, SWC, G,...
                        GPP, ET, Site, Type, ....
                        MAH_all, obs_length,...
                        WUE);
            elseif (mah==2)
                fig_desp=['_MAHbase'];
                X = table(Tair, VPD, CO2, SW, USTAR, Tsoil, SWC, G,...
                        GPP, ET, Site, Type, ....
                        MAH_base, obs_length,...
                        WUE);
            elseif (mah==3)
                fig_desp=['_MAHtran'];
                X = table(Tair, VPD, CO2, SW, USTAR, Tsoil, SWC, G,...
                        GPP, ET, Site, Type, ....
                        MAH_tran, obs_length,...
                        WUE);
            end
            t = templateTree('NumPredictorsToSample','all', ...
                'PredictorSelection', 'curvature', 'Surrogate', 'on');
            rng('default'); % For reproducibility
            Mdl = fitrensemble(X, 'WUE', 'Method', 'bag', 'NumLearningCycles', 100, 'Learners', t);
            imp = oobPermutedPredictorImportance(Mdl);
            subplot(2,2,variable)
            [tmp, idx] = sort(imp,'descend');
            col = length(tmp);
            bar(tmp(1:col));
            xlabel('Predictors', 'FontSize', 12);
            % adjust xtick format
            var_name = Mdl.PredictorNames;
            var_name{5} = '{\it u*}';
            if (mah==1)
                var_name{13} = '{\it MAH_{{78-99}}}';
            elseif (mah==2)
                var_name{13} = '{\it MAH_{{40-60}}}';
            elseif (mah==3)
                var_name{13} = '{\it MAH_{{51-77}}}';
            end
            var_name{14} = '{Obs length}';
            set(gca,'xtick',[1:col])
            set(gca,'xticklabel',string(var_name(idx)), 'XTickLabelRotation', -45)
            ylim([0 max(imp)+1])
            ylabel('Predictor importance estimates', 'FontSize', 12);
            xlabel('Normalized IAV (% yr^{-1})', 'FontSize', 12);
            id = char(erase(id_name(variable),'"'));
            t = text(0.02,0.98,[id ' ' yname ' IAV'],...
                'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12);
            subplot(2,2,variable+2)
            % R2
            yHat = oobPredict(Mdl);
            R2 = corr(Mdl.Y,yHat)^2;
            % one to one
            plot([min(Mdl.Y),max(Mdl.Y)],[min(Mdl.Y),max(Mdl.Y)],'k--')
            hold on
            % line fit
            rg_fit_tmp = polyfit(Mdl.Y, yHat, 1);
            xxx = linspace(min(Mdl.Y),max(Mdl.Y),10);
            fff = polyval(rg_fit_tmp,xxx);
            plot(xxx,fff,'b-','Linewidth',2)
            ppp = dscatter(Mdl.Y, yHat);
            colormap('Hot')
            axis([min(Mdl.Y) max(Mdl.Y) min(Mdl.Y) max(Mdl.Y)])
            xlabel(['Observed ' yname ' IAV (% yr^{-1})'], 'FontSize', 12);
            ylabel(['Estimated ' yname ' IAV (% yr^{-1})'], 'FontSize', 12);
            id = char(erase(id_name(variable+2),'"'));
            t = text(0.02,0.98,[id ' ' yname],...
                'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                    'FontSize',12);
            t = text(0.02,0.78,['R^2 = ' num2str(sprintf('%1.2f',R2))],...
                'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                    'FontSize',12);
        end
        figname=['../plots_hourly/uWUE_trend_full_v7/' fig_header  fig_desp '.jpeg'];
        print('-djpeg','-r300',[figname]);
    end    
end

%% what makes the difference between MAH_all and MAH_base?
% sites with high MAH_all importance and high MAH_base importance
idx_high = find(r_Sen_MAH_pct_base_uwue>=0.5);
% sites with high MAH_all importance and low MAH_base importance
idx_low = find(r_Sen_MAH_pct_base_uwue<0.5);
clc;
close all
units = '(g C hPa^{0.5} kg H_2O^{-1})';
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", ...
    "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", ...
    "(q)", "(r)", "(s)", "(t)", "(u)", "(v)", "(w)", "(x)"];

for r_val = 1:2
    clf;
    if (r_val==1)
        idx_site = idx_high;
        fig_type = '_r_high';
    elseif (r_val==2)
        idx_site = idx_low;
        fig_type = '_r_low';
    end
    for type = 1:1
        clf;
        if (type==1)
            wue_data = uwue_MAH_grs_mean_subset;
            wue_name = '{\it uWUE}';
            wue_med_lgd = '{\it uWUE_{median}}';
            wue_agg_lgd = '{\it uWUE_{aggregate}}';
            fig_header = 'RF_VR_WUE_IAV_uWUE_grs';
        elseif (type==2)
            wue_data = iwue_MAH_grs_mean_subset;
            wue_name = '{\it WUE_{ei}}';
            wue_med_lgd = '{\it WUE_{ei, median}}';
            wue_agg_lgd = '{\it WUE_{ei, aggregate}}';
            fig_header = 'RF_VR_WUE_IAV_iWUE_grs';
        elseif (type==3)
            wue_data = uwue_MAH_sum_mean_subset;
            wue_name = '{\it uWUE}';
            wue_med_lgd = '{\it uWUE_{median}}';
            wue_agg_lgd = '{\it uWUE_{aggregate}}';
            fig_header = 'RF_VR_WUE_IAV_uWUE_sum';
        elseif (type==4)
            wue_data = iwue_MAH_sum_mean_subset;
            wue_name = '{\it WUE_{ei}}';
            wue_med_lgd = '{\it WUE_{ei, median}}';
            wue_agg_lgd = '{\it WUE_{ei, aggregate}}';
            fig_header = 'RF_VR_WUE_IAV_iWUE_sum';
        end
        % normalized trend estimates
        grs_trend_agg = [];
        for i = 1:length(idx_site)
            idx = find(wue_data(:, 19)==site_MAH_ana(idx_site(i)));
            data_site = wue_data(idx,:);
            tmp = (data_site(2:end,:)-data_site(1:end-1,:))./...
                (data_site(2:end,17)-data_site(1:end-1,17))./...
                mean(data_site,'omitnan');
            % include site and type info
            if (length(tmp)>0)
                tmp(:,19) = site_included(i);
                tmp(:,39) = data_site(1, 39);
            end
            grs_trend_agg = [grs_trend_agg; tmp];
        end
        % convert to % per year
        grs_trend_agg = grs_trend_agg*100;
        for variable = 1:2
            if (variable==1)
                yname = wue_med_lgd;
                WUE = grs_trend_agg(:,1);
            elseif (variable==2)
                yname = wue_agg_lgd;
                WUE = grs_trend_agg(:,34);
            end
            idx = find(isfinite(WUE)==1);
            WUE = WUE(idx);
            Tair = grs_trend_agg(idx,3);
            VPD = grs_trend_agg(idx,5);
            CO2 = grs_trend_agg(idx,6);
            SW = grs_trend_agg(idx,7);
            USTAR = grs_trend_agg(idx,8);
            Wind = grs_trend_agg(idx,9);
            PPFD = grs_trend_agg(idx,10);
            Tsoil = grs_trend_agg(idx,11);
            SWC = grs_trend_agg(idx,12);
            G = grs_trend_agg(idx,13);
            NEE = grs_trend_agg(idx,14);
            GPP = grs_trend_agg(idx,15);
            ET = grs_trend_agg(idx,16);
            MAH_all = grs_trend_agg(idx,35);
            MAH_base = grs_trend_agg(idx,36);
            obs_length = grs_trend_agg(idx,38);
            DOY = categorical(round(grs_trend_agg(idx,2)));
            Hour = categorical(round(grs_trend_agg(idx,18)));
            Site = categorical(grs_trend_agg(idx,19));
            Type = categorical(grs_trend_agg(idx,39));
            X = table(Tair, VPD, CO2, SW, USTAR, Tsoil, SWC, G,...
                GPP, ET, Site, Type, ....
                MAH_all, MAH_base, obs_length,...
                WUE);
            t = templateTree('NumPredictorsToSample','all', ...
                'PredictorSelection', 'curvature', 'Surrogate', 'on');
            rng('default'); % For reproducibility
            Mdl = fitrensemble(X, 'WUE', 'Method', 'bag', 'NumLearningCycles', 200, 'Learners', t);
            imp = oobPermutedPredictorImportance(Mdl);
            subplot(2,2,variable)
            [tmp, idx] = sort(imp,'descend');
            col = length(tmp);
            bar(tmp(1:col));
            xlabel('Predictors', 'FontSize', 12);
            % adjust xtick format
            var_name = Mdl.PredictorNames;
            var_name{5} = '{\it u*}';
            var_name{13} = '{\it MAH_{{78-99}}}';
            var_name{14} = '{\it MAH_{{40-60}}}';
            var_name{15} = '{Obs length}';
            set(gca,'xtick',[1:col])
            set(gca,'xticklabel',string(var_name(idx)), 'XTickLabelRotation', -45)
            ylim([0 max(imp)+1])
            ylabel('Predictor importance estimates', 'FontSize', 12);
            xlabel('Normalized IAV (% yr^{-1})', 'FontSize', 12);
            id = char(erase(id_name(variable),'"'));
            t = text(0.02,0.98,[id ' ' yname ' IAV'],...
                'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12);
            subplot(2,2,variable+2)
            % R2
            yHat = oobPredict(Mdl);
            R2 = corr(Mdl.Y,yHat)^2;
            % one to one
            plot([min(Mdl.Y),max(Mdl.Y)],[min(Mdl.Y),max(Mdl.Y)],'k--')
            hold on
            % line fit
            rg_fit_tmp = polyfit(Mdl.Y, yHat, 1);
            xxx = linspace(min(Mdl.Y),max(Mdl.Y),10);
            fff = polyval(rg_fit_tmp,xxx);
            plot(xxx,fff,'b-','Linewidth',2)
            ppp = dscatter(Mdl.Y, yHat);
            colormap('Hot')
        %     c = colorbar('Position', [0.92, 0.17, 0.02, 0.3]);
            axis([min(Mdl.Y) max(Mdl.Y) min(Mdl.Y) max(Mdl.Y)])
            xlabel(['Observed ' yname ' IAV (% yr^{-1})'], 'FontSize', 12);
            ylabel(['Estimated ' yname ' IAV (% yr^{-1})'], 'FontSize', 12);
            id = char(erase(id_name(variable+2),'"'));
            t = text(0.02,0.98,[id ' ' yname],...
                'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                    'FontSize',12);
            t = text(0.02,0.78,['R^2 = ' num2str(sprintf('%1.2f',R2))],...
                'Units', 'Normalized', 'VerticalAlignment', 'Top',...
                    'FontSize',12);
        end
        figname=['../plots_hourly/uWUE_trend_full_v7/' fig_header fig_type '.jpeg'];
        print('-djpeg','-r300',[figname]);
    end
end


%% examine factors controlling MAH, based on climates labeled at individual WUE percentiles
% only include the measued data
clc;
close all;
% locate the variables: Tair; Precip; VPD; CO2 mole fraction; GPP; 
% Solar radiation; USTAR
var_test = ["TA_F","P_F","VPD_F",...
    "CO2_F_MDS","SW_IN_F_MDS","USTAR","WS_F","PPFD_IN","TS_F_MDS_1",...
    "SWC_F_MDS_1", "G_F_MDS","NEE_VUT_REF"];
var_idx = [];
for i = 1:length(var_test)
    var_idx = [var_idx; find(strcmpi(var_test(i),var_wanted))];
end
% measured data only
data = [];
for j = 1:length(var_idx)
    col = var_idx(j);
    data(:,j+1) = hourly_data(:,col);
    % QC flags, only use measured data
    idx = find(hourly_data(:,col+1)>0);
    data(idx,j+1) = nan;
end
% set time label
data(:,1) = Time(:,2);
% locate GPP
% umolCO2 m-2 s-1 for hourly data
var_test = ["GPP_NT_VUT_REF"];
var_idx = [];
for i = 1:length(var_test)
    var_idx = [var_idx; find(strcmpi(var_test(i),var_wanted))];
end
% GPP QC
gpp = hourly_data(:,var_idx(1));
% set GPP, umolCO2 m-2 s-1 to gCm-2s-1
data(:,end+1) = gpp*12*(10^-6);
% set ET
data(:,end+1) = daily_ET;
% include site-year info
data(:,end+1) = Time(:,1);
% include site-hour info
data(:,end+1) = Time(:,3);
% include site info
data(:,end+1) = hourly_data(:,end);
vpd_m = [];
Tair_m = [];
cco2_m = [];
% define daytime as solar radiation > 100 W/m2
sw_threshold = 100;
et_threshold = 10^-6;
gpp_threshold = 10^-6;
tair_threshold = 0;
ustar_threshold = 0.1;
gpp_pct = [];
et_pct = [];
vpd_pct = [];
uwue_MAH_control_grs = [];
iwue_MAH_control_grs = [];
% declare month DOY
for month = 2:12
    month_start(month) = mon_sum(month-1)+1;
    month_end(month) = mon_sum(month);
end
month_start(1) = 1;
month_end(1) = mon_sum(1);
pct_window = 1:99;
pct_window_2pct = 1:2:99;
clear uwue uwue obs_year pct_value uwue_threshold_bins iwue_threshold_bins
% calculate WUE based on seasonal means/sums (1 WUE per growing season)
for site = 1:length(site_included)
    % exclude croplands
    if (IGBP(site)~='CRO'|IGBP(site)~='WET')
        count = 1;
        idx = find(data(:,18)==site_included(site));
        data_site = data(idx,:);
        yr_range = [data_site(1, 16):data_site(end, 16)];
        if (length(yr_range)>0)
            yr_start(site) = yr_range(1);
            yr_end(site) = yr_range(end);
            % WUE for the entire time series
            idx = find(wue_hourly_agg(:,4)==site_included(site));
            if (length(idx)>0)
                % uwue
                tmp = wue_hourly_agg(idx, 3);
                uwue_threshold = [prctile(tmp,pct_best_r_agg(1)),prctile(tmp,pct_best_r_agg(1)+1)];
                % 1 percent windows
                for bin = 1:length(pct_window)
                    uwue_threshold_bins(bin) = prctile(tmp,pct_window(bin));
                end
                % iwue
                tmp = wue_hourly_agg(idx, 15);
                iwue_threshold = [prctile(tmp,pct_best_r_agg(2)),prctile(tmp,pct_best_r_agg(2)+1)];
                % 1 percent windows
                for bin = 1:length(pct_window)
                    iwue_threshold_bins(bin) = prctile(tmp,pct_window(bin));
                end
            end
            for yr = 1:length(yr_range)
                idx = find(data_site(:,16)==yr_range(yr));
                data_site_year = data_site(idx,:);
                % precip events
                precip = sum(data_site_year(:,3),'omitnan'); % total precip [mm]
                % precip events
                idx = find(data_site_year(:,3)>0);
                doy_tmp = unique(data_site_year(idx,1));
                % filter out data during and one day after a precip event
                for dd = 1:length(doy_tmp)
                    idx = find(data_site_year(:,1)==doy_tmp(dd)|...
                        data_site_year(:,1)==(doy_tmp(dd)+1));
                    data_site_year(idx,:) = nan;
                end
                % filter out low ustar
                idx = find(data_site_year(:,7)<ustar_threshold);
                data_site_year(idx,:) = nan;
                % filter out nighttime data
                idx = find(data_site_year(:,6)<sw_threshold);
                data_site_year(idx,:) = nan;
                % condensation events
                idx = find(data_site_year(:,15)<et_threshold);
                data_site_year(idx,:) = nan;
                % set growing season as the period when daily GPP>10%
                % daily GPP max and daily Tair>0
                idx = isfinite(data_site_year(:, 1));
                days = unique(data_site_year(idx, 1));
                tmp = [];
                for day = 1:length(days)
                    idx = find(data_site_year(:, 1)==days(day));
                    if (length(idx)>0)
                        % DOY, Tair, GPP
                        xxx = [data_site_year(idx(1), 1), ...
                            mean(data_site_year(idx,2), 'omitnan'),...
                            sum(data_site_year(idx,14), 'omitnan')];
                        tmp = [tmp; xxx];
                    end
                end
                idx = find(tmp(:,2)<=0 | tmp(:,3)<max(tmp(:,3))*0.05);
                if (length(idx)>1)
                    days = unique(tmp(idx, 1));
                    for day = 1:length(days)
                        idx = find(data_site_year(:, 1)==days(day));
                        if (length(idx)>0)
                            data_site_year(idx,:) = nan;
                        end
                    end
                end
                % When GPP<0
                idx = find(data_site_year(:,14)<gpp_threshold);
                data_site_year(idx,:) = nan;
                % When Tair<0
                idx = find(data_site_year(:,2)<tair_threshold);
                data_site_year(idx,:) = nan;
                uwue_hourly = data_site_year(:,14).*(data_site_year(:,4).^.5)./data_site_year(:,15);
                iwue_hourly = data_site_year(:,14).*(data_site_year(:,4))./data_site_year(:,15);
                % daily aggregate
                % output year, DOY, WUE, site, GPP, ET, VPD, Tair, [CO2], USTAR
                idx = isfinite(data_site_year(:, 1));
                days = unique(data_site_year(idx, 1));
                tmp = [];
                for day = 1:length(days)
                    idx = find(data_site_year(:, 1)==days(day));
                    if (length(idx)>0)
                        xxx = [data_site_year(idx(1),10), data_site_year(idx(1), 1), ...
                        mean(uwue_hourly(idx), 'omitnan'), data_site_year(idx(1),12), ...
                        sum(data_site_year(idx,14), 'omitnan'), sum(data_site_year(idx,15), 'omitnan'), ...
                        mean(data_site_year(idx,4), 'omitnan'), mean(data_site_year(idx,2), 'omitnan'), ...
                        mean(data_site_year(idx,5), 'omitnan'), mean(data_site_year(idx,7), 'omitnan')];
                        xxx(3) = xxx(5)*(xxx(7)^.5)/xxx(6); % WUE from daily agg values
                        tmp = [tmp; xxx];
                    end
                end
                uwue_daily = tmp;
                if (length(uwue_daily)>0) % data exist
                    if (sum(isfinite(uwue_daily(:,3)))>=1)
                        % uWUE
                        tmp = [uwue_hourly, data_site_year];
                        % categorize by WUE bins
                        for bin = 1:length(uwue_threshold_bins)
                            if (bin==length(uwue_threshold_bins))
                                idx = find(uwue_hourly>=uwue_threshold_bins(bin));
                            else
                                idx = find(uwue_hourly>=uwue_threshold_bins(bin) & uwue_hourly<uwue_threshold_bins(bin+1));
                            end
                            if (length(idx)>0)
                                output = [tmp(idx,:), ones(length(idx),1).*pct_window(bin)];
                            else
                                output = [];
                            end
                            uwue_MAH_control_grs = [uwue_MAH_control_grs; output];
                        end
                        % WUEei
                        tmp = [iwue_hourly, data_site_year];
                        % categorize by WUE bins
                        for bin = 1:length(iwue_threshold_bins)
                            if (bin==length(iwue_threshold_bins))
                                idx = find(iwue_hourly>=iwue_threshold_bins(bin));
                            else
                                idx = find(iwue_hourly>=iwue_threshold_bins(bin) & iwue_hourly<iwue_threshold_bins(bin+1));
                            end
                            if (length(idx)>0)
                                output = [tmp(idx,:), ones(length(idx),1).*pct_window(bin)];
                            else
                                output = [];
                            end
                            iwue_MAH_control_grs = [iwue_MAH_control_grs; output];
                        end
                    end
                end
            end
        end
    end
end
% filterout nan values
uwue_MAH_control_grs(uwue_MAH_control_grs==0) = NaN;
iwue_MAH_control_grs(iwue_MAH_control_grs==0) = NaN;

%% examine factors controlling MAH, based on climates labeled at individual WUE percentiles
% only include the measued data
clc;
close all;
% locate the variables: Tair; Precip; VPD; CO2 mole fraction; GPP; 
% Solar radiation; USTAR
var_test = ["TA_F","P_F","VPD_F",...
    "CO2_F_MDS","SW_IN_F_MDS","USTAR","WS_F","PPFD_IN","TS_F_MDS_1",...
    "SWC_F_MDS_1", "G_F_MDS","NEE_VUT_REF"];
var_idx = [];
for i = 1:length(var_test)
    var_idx = [var_idx; find(strcmpi(var_test(i),var_wanted))];
end
% measured data only
data = [];
for j = 1:length(var_idx)
    col = var_idx(j);
    data(:,j+1) = hourly_data(:,col);
    % QC flags, only use measured data
    idx = find(hourly_data(:,col+1)>0);
    data(idx,j+1) = nan;
end
% set time label
data(:,1) = Time(:,2);
% locate GPP
% umolCO2 m-2 s-1 for hourly data
var_test = ["GPP_NT_VUT_REF"];
var_idx = [];
for i = 1:length(var_test)
    var_idx = [var_idx; find(strcmpi(var_test(i),var_wanted))];
end
% GPP QC
gpp = hourly_data(:,var_idx(1));
% set GPP, umolCO2 m-2 s-1 to gCm-2s-1
data(:,end+1) = gpp*12*(10^-6);
% set ET
data(:,end+1) = daily_ET;
% include site-year info
data(:,end+1) = Time(:,1);
% include site-hour info
data(:,end+1) = Time(:,3);
% include site info
data(:,end+1) = hourly_data(:,end);
vpd_m = [];
Tair_m = [];
cco2_m = [];
% define daytime as solar radiation > 100 W/m2
sw_threshold = 100;
et_threshold = 10^-6;
gpp_threshold = 10^-6;
tair_threshold = 0;
ustar_threshold = 0.1;
uwue_MAH_control_sum = [];
iwue_MAH_control_sum = [];
% declare month DOY
for month = 2:12
    month_start(month) = mon_sum(month-1)+1;
    month_end(month) = mon_sum(month);
end
month_start(1) = 1;
month_end(1) = mon_sum(1);
pct_window = 1:99;
pct_window_2pct = 1:2:99;
clear uwue uwue obs_year pct_value uwue_threshold_bins iwue_threshold_bins
% calculate WUE based on seasonal means/sums (1 WUE per growing season)
for site = 1:length(site_included)
    lat_tmp = LAT(site);
    % exclude croplands
    if (IGBP(site)~='CRO'|IGBP(site)~='WET')
        count = 1;
        idx = find(data(:,18)==site_included(site));
        data_site = data(idx,:);
        yr_range = [data_site(1, 16):data_site(end, 16)];
        if (length(yr_range)>0)
            yr_start(site) = yr_range(1);
            yr_end(site) = yr_range(end);
            % WUE for the entire time series
            idx = find(wue_hourly_agg(:,4)==site_included(site));
            if (length(idx)>0)
                % uwue
                tmp = wue_hourly_agg(idx, 3);
                uwue_threshold = [prctile(tmp,pct_best_r_agg(1)),prctile(tmp,pct_best_r_agg(1)+1)];
                % 1 percent windows
                for bin = 1:length(pct_window)
                    uwue_threshold_bins(bin) = prctile(tmp,pct_window(bin));
                end
                % iwue
                tmp = wue_hourly_agg(idx, 15);
                iwue_threshold = [prctile(tmp,pct_best_r_agg(2)),prctile(tmp,pct_best_r_agg(2)+1)];
                % 1 percent windows
                for bin = 1:length(pct_window)
                    iwue_threshold_bins(bin) = prctile(tmp,pct_window(bin));
                end
            end
            for yr = 1:length(yr_range)
                idx = find(data_site(:,16)==yr_range(yr));
                data_site_year = data_site(idx,:);
                % precip events
                precip = sum(data_site_year(:,3),'omitnan'); % total precip [mm]
                % precip events
                idx = find(data_site_year(:,3)>0);
                doy_tmp = unique(data_site_year(idx,1));
                % filter out data during and one day after a precip event
                for dd = 1:length(doy_tmp)
                    idx = find(data_site_year(:,1)==doy_tmp(dd)|...
                        data_site_year(:,1)==(doy_tmp(dd)+1));
                    data_site_year(idx,:) = nan;
                end
                % filter out low ustar
                idx = find(data_site_year(:,7)<ustar_threshold);
                data_site_year(idx,:) = nan;
                % filter out nighttime data
                idx = find(data_site_year(:,6)<sw_threshold);
                data_site_year(idx,:) = nan;
                % condensation events
                idx = find(data_site_year(:,15)<et_threshold);
                data_site_year(idx,:) = nan;
                % define summer season
                if (lat_tmp>0)
                    summer_start = mon_sum(5)+1;
                    summer_end = mon_sum(8);
                    idx_other = find(data_site_year(:, 1)<summer_start | ...
                        data_site_year(:, 1)>summer_end);
                else
                    summer_start = mon_sum(11)+1;
                    summer_end = mon_sum(2);
                    idx_other = find(data_site_year(:, 1)<summer_start & ...
                        data_site_year(:, 1)>summer_end);
                end
                data_site_year(idx_other,:) = nan;
                idx = isfinite(data_site_year(:, 1));
                days = unique(data_site_year(idx, 1));
                tmp = [];
                for day = 1:length(days)
                    idx = find(data_site_year(:, 1)==days(day));
                    if (length(idx)>0)
                        % DOY, Tair, GPP
                        xxx = [data_site_year(idx(1), 1), ...
                            mean(data_site_year(idx,2), 'omitnan'),...
                            sum(data_site_year(idx,14), 'omitnan')];
                        tmp = [tmp; xxx];
                    end
                end
                idx = find(tmp(:,2)<=0 | tmp(:,3)<max(tmp(:,3))*0.05);
                if (length(idx)>1)
                    days = unique(tmp(idx, 1));
                    for day = 1:length(days)
                        idx = find(data_site_year(:, 1)==days(day));
                        if (length(idx)>0)
                            data_site_year(idx,:) = nan;
                        end
                    end
                end
                % When GPP<0
                idx = find(data_site_year(:,14)<gpp_threshold);
                data_site_year(idx,:) = nan;
                % When Tair<0
                idx = find(data_site_year(:,2)<tair_threshold);
                data_site_year(idx,:) = nan;
                uwue_hourly = data_site_year(:,14).*(data_site_year(:,4).^.5)./data_site_year(:,15);
                iwue_hourly = data_site_year(:,14).*(data_site_year(:,4))./data_site_year(:,15);
                % daily aggregate
                % output year, DOY, WUE, site, GPP, ET, VPD, Tair, [CO2], USTAR
                idx = isfinite(data_site_year(:, 1));
                days = unique(data_site_year(idx, 1));
                tmp = [];
                for day = 1:length(days)
                    idx = find(data_site_year(:, 1)==days(day));
                    if (length(idx)>0)
                        xxx = [data_site_year(idx(1),10), data_site_year(idx(1), 1), ...
                        mean(uwue_hourly(idx), 'omitnan'), data_site_year(idx(1),12), ...
                        sum(data_site_year(idx,14), 'omitnan'), sum(data_site_year(idx,15), 'omitnan'), ...
                        mean(data_site_year(idx,4), 'omitnan'), mean(data_site_year(idx,2), 'omitnan'), ...
                        mean(data_site_year(idx,5), 'omitnan'), mean(data_site_year(idx,7), 'omitnan')];
                        xxx(3) = xxx(5)*(xxx(7)^.5)/xxx(6); % WUE from daily agg values
                        tmp = [tmp; xxx];
                    end
                end
                uwue_daily = tmp;
                if (length(uwue_daily)>0) % data exist
                    % output year, DOY, WUE, site, GPP, ET, VPD, Tair, [CO2], USTAR
                    if (sum(isfinite(uwue_daily(:,3)))>=1)
                        % uWUE
                        tmp = [uwue_hourly, data_site_year];
                        % categorize by WUE bins
                        for bin = 1:length(uwue_threshold_bins)
                            if (bin==length(uwue_threshold_bins))
                                idx = find(uwue_hourly>=uwue_threshold_bins(bin));
                            else
                                idx = find(uwue_hourly>=uwue_threshold_bins(bin) & uwue_hourly<uwue_threshold_bins(bin+1));
                            end
                            if (length(idx)>0)
                                output = [tmp(idx,:), ones(length(idx),1).*pct_window(bin)];
                            else
                                output = [];
                            end
                            uwue_MAH_control_sum = [uwue_MAH_control_sum; output];
                        end
                        % WUEei
                        tmp = [iwue_hourly, data_site_year];
                        % categorize by WUE bins
                        for bin = 1:length(iwue_threshold_bins)
                            if (bin==length(iwue_threshold_bins))
                                idx = find(iwue_hourly>=iwue_threshold_bins(bin));
                            else
                                idx = find(iwue_hourly>=iwue_threshold_bins(bin) & iwue_hourly<iwue_threshold_bins(bin+1));
                            end
                            if (length(idx)>0)
                                output = [tmp(idx,:), ones(length(idx),1).*pct_window(bin)];
                            else
                                output = [];
                            end
                            iwue_MAH_control_sum = [iwue_MAH_control_sum; output];
                        end
                    end
                end
            end
        end
    end
end
% filterout nan values
uwue_MAH_control_sum(uwue_MAH_control_sum==0) = NaN;
iwue_MAH_control_sum(iwue_MAH_control_sum==0) = NaN;

%% analyze results when there are moren than 10 years of measurements
clc;
IGBP_included = unique(IGBP);% % exclude croplands
idx = find(IGBP_included=='CRO'|IGBP_included=='WET');
IGBP_included(idx) = [];
uwue_MAH_pct_ctrl_grs_subset = [];
iwue_MAH_pct_ctrl_grs_subset = [];
uwue_MAH_pct_ctrl_sum_subset = [];
iwue_MAH_pct_ctrl_sum_subset = [];
% uwue_MAH_pct_subset = [];
for i = 1:length(IGBP_included)
    % find sites within an ecosystem type
    idx_IGBP = find(IGBP==IGBP_included(i));
    for j = 1:length(idx_IGBP)
        site = site_included(idx_IGBP(j));
        % extract site data
        % uwue, grs
        idx = find(uwue_MAH_control_grs(:,19)==site);
        site_year = unique(uwue_MAH_control_grs(idx,17));
        % include the site data if there's more than 10 site-years
        if (length(site_year)>=10)
            tmp = [uwue_MAH_control_grs(idx,:), ones(length(idx),1)*i]; % attach type info
            uwue_MAH_pct_ctrl_grs_subset = [uwue_MAH_pct_ctrl_grs_subset; tmp];
        end
        % iwue, grs
        idx = find(iwue_MAH_control_grs(:,19)==site);
        site_year = unique(iwue_MAH_control_grs(idx,17));
        % include the site data if there's more than 10 site-years
        if (length(site_year)>=10)
            tmp = [iwue_MAH_control_grs(idx,:), ones(length(idx),1)*i]; % attach type info
            iwue_MAH_pct_ctrl_grs_subset = [iwue_MAH_pct_ctrl_grs_subset; tmp];
        end
        % uwue, summer
        idx = find(uwue_MAH_control_sum(:,19)==site);
        site_year = unique(uwue_MAH_control_sum(idx,17));
        % include the site data if there's more than 10 site-years
        if (length(site_year)>=10)
            tmp = [uwue_MAH_control_sum(idx,:), ones(length(idx),1)*i]; % attach type info
            uwue_MAH_pct_ctrl_sum_subset = [uwue_MAH_pct_ctrl_sum_subset; tmp];
        end
        % iwue, summer
        idx = find(iwue_MAH_control_sum(:,19)==site);
        site_year = unique(iwue_MAH_control_sum(idx,17));
        % include the site data if there's more than 10 site-years
        if (length(site_year)>=10)
            tmp = [iwue_MAH_control_sum(idx,:), ones(length(idx),1)*i]; % attach type info
            iwue_MAH_pct_ctrl_sum_subset = [iwue_MAH_pct_ctrl_sum_subset; tmp];
        end
    end
end
%% calculate annual changes in individual uWUE percentiles
clc;
sites = unique(uwue_MAH_pct_ctrl_grs_subset(:,19));
pct_window = 1:99;
uwue_MAH_pct_ctrl_median_grs = [];
iwue_MAH_pct_ctrl_median_grs = [];
uwue_MAH_pct_ctrl_median_sum = [];
iwue_MAH_pct_ctrl_median_sum = [];
for type = 1:4
    if (type==1)
        wue_data = uwue_MAH_pct_ctrl_grs_subset;
    elseif (type==2)
        wue_data = iwue_MAH_pct_ctrl_grs_subset;
    elseif (type==3)
        wue_data = uwue_MAH_pct_ctrl_sum_subset;
    elseif (type==4)
        wue_data = iwue_MAH_pct_ctrl_sum_subset;
    end
    for i = 1:length(sites)
        pct_median_site = [];
        idx = find(wue_data(:,19)==sites(i));
        data_site = wue_data(idx,:);
        years = unique(data_site(:,17));
        for j = 1:length(years)
            idx = find(data_site(:,17)==years(j));
            data_site_year = data_site(idx,:);
            for k = 1:length(pct_window)
                idx = find(data_site_year(:,20)==pct_window(k));
                if (length(idx)>1)
                    tmp_median = median(data_site_year(idx,:), 'omitnan');
                else
                    tmp_median = data_site_year(idx,:);
                end
                pct_median_site = [pct_median_site; tmp_median];
            end
        end
        if (type==1)
            uwue_MAH_pct_ctrl_median_grs = [uwue_MAH_pct_ctrl_median_grs; pct_median_site];
        elseif (type==2)
            iwue_MAH_pct_ctrl_median_grs = [iwue_MAH_pct_ctrl_median_grs; pct_median_site];
        elseif (type==3)
            uwue_MAH_pct_ctrl_median_sum = [uwue_MAH_pct_ctrl_median_sum; pct_median_site];
        elseif (type==4)
            iwue_MAH_pct_ctrl_median_sum = [iwue_MAH_pct_ctrl_median_sum; pct_median_site];
        end        
    end
end

%% random-forest model variable importance: which factor controls WUE at individual percentiles?
% are the factors contributing to uWUE sensitive to ecosystem types?
close all;
clc;
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", ...
    "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", ...
    "(q)", "(r)", "(s)", "(t)", "(u)", "(v)", "(w)", "(x)"];
for type = 1:1
    if (type==1)
        wue_data = uwue_MAH_pct_ctrl_median_grs;
        path = '../plots_hourly/uWUE_trend_full_v7/uwue_MAH_pct_ctrl_grs/';
        wue_name = '{\it uWUE}';
        sup1 = 'th';
        sup2 = 'th';
    elseif (type==2)
        wue_data = iwue_MAH_pct_ctrl_median_grs;
        path = '../plots_hourly/uWUE_trend_full_v7/iwue_MAH_pct_ctrl_grs/';
        wue_name = '{\it WUE_{ei}}';
        sup1 = 'th';
        sup2 = 'th';
    elseif (type==3)
        wue_data = uwue_MAH_pct_ctrl_median_sum;
        path = '../plots_hourly/uWUE_trend_full_v7/uwue_MAH_pct_ctrl_sum/';
        wue_name = '{\it uWUE}';
        sup1 = 'th';
        sup2 = 'th';
    elseif (type==4)
        wue_data = iwue_MAH_pct_ctrl_median_sum;
        path = '../plots_hourly/uWUE_trend_full_v7/iwue_MAH_pct_ctrl_sum/';
        wue_name = '{\it WUE_{ei}}';
        sup1 = 'rd';
        sup2 = 'nd';
    end
    % group into forest and non-forest sites
    for pft = 1:3
        tmp = wue_data(:,21);
        if (pft==1) % forest
            idx_IGBP = find(tmp==2|tmp==3|tmp==4|tmp==6);
            fig_type = '_forest';
        elseif (pft==2)
            idx_IGBP = find(tmp==5|tmp==8);
            fig_type = '_n_forest';
        elseif (pft==3) % all
            idx_IGBP = find(tmp==2|tmp==3|tmp==4|tmp==5|tmp==6|tmp==8);
            fig_type = '_all';
        end
        r2_lower_pct = [];
        vr_lower_pct = [];
        dp_lower_pct = [];
        r2_higher_pct = [];
        vr_higher_pct = [];
        dp_higher_pct = [];
        r2_all_pct = [];
        vr_all_pct = [];
        dp_all_pct = [];
        site_vr_agg_all = [];
        site_vr_agg_high = [];
        site_vr_agg_low = [];
        sites = unique(wue_data(idx_IGBP,19));
        for site = 1:length(sites)
            clf;
            siteID = dir(sites(site), 5:10);
            for variable = 1:3
                if (variable==1)
                    idx = find(wue_data(:,19)==sites(site) & wue_data(:,20)<=99);
                    desp = [wue_name '_{1^{st}-99^{th}}'];
                elseif (variable==2)
                    idx = find(wue_data(:,19)==sites(site) & wue_data(:,20)>=pct_best_r_agg(type) & wue_data(:,20)<=99);
                    desp = [wue_name '_{' int2str(pct_best_r_agg(type)) '^{' sup1 '}-99^{th}}'];
                elseif (variable==3)
                    idx = find(wue_data(:,19)==sites(site) & wue_data(:,20)>=1 & wue_data(:,20)<pct_best_r_agg(type));
                    desp = [wue_name '_{1^{st}-' int2str(pct_best_r_agg(type)-1) '^{' sup2 '}}'];
                end
                % build up variable list
                nt = length(idx);
                Tair = wue_data(idx,3);
                VPD = wue_data(idx,5);
                CO2 = wue_data(idx,6);
                SW = wue_data(idx,7);
                USTAR = wue_data(idx,8);
                Wind = wue_data(idx,9);
                PPFD = wue_data(idx,10);
                Tsoil = wue_data(idx,11);
                SWC = wue_data(idx,12);
                G = wue_data(idx,13);
                GPP = wue_data(idx,15);
                ET = wue_data(idx,16);
                WUE = wue_data(idx,1);
                X = table(Tair, VPD, CO2, SW, USTAR, Tsoil, SWC, G, ...
                    WUE);
                t = templateTree('NumPredictorsToSample','all', ...
                    'PredictorSelection', 'curvature', 'Surrogate', 'on');
                rng('default'); % For reproducibility
                Mdl = fitrensemble(X, 'WUE', 'Method', 'bag', 'NumLearningCycles', 100, 'Learners', t);
                imp = oobPermutedPredictorImportance(Mdl);
                [tmp, idx] = sort(imp,'descend');

                % store output for later comparison
                if (variable==3)
                    vr_lower_pct = [vr_lower_pct; imp];
                    r2_lower_pct = [r2_lower_pct; R2];
                    dp_lower_pct = [dp_lower_pct; nt];
                    site_vr_agg_low = [site_vr_agg_low; idx];
                elseif (variable==2)
                    vr_higher_pct = [vr_higher_pct; imp];
                    r2_higher_pct = [r2_higher_pct; R2];
                    dp_higher_pct = [dp_higher_pct; nt];
                    site_vr_agg_high = [site_vr_agg_high; idx];
                elseif (variable==1)
                    vr_all_pct = [vr_all_pct; imp];
                    r2_all_pct = [r2_all_pct; R2];
                    dp_all_pct = [dp_all_pct; nt];
                    site_vr_agg_all = [site_vr_agg_all; idx];
                end
                vr_name = Mdl.PredictorNames;
                vr_name{5} = 'u*';
            end
        end
        % RF variable importance across sites
        clc;
        close all
        id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", ...
            "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", ...
            "(q)", "(r)", "(s)", "(t)", "(u)", "(v)", "(w)", "(x)"];

        for variable = 1:3
            if (variable==1)
                tmp_mean = mean(vr_all_pct);
                tmp_std = std(vr_all_pct);
                nt = [mean(dp_all_pct), std(dp_all_pct), sum(dp_all_pct)];
                desp = [wue_name '_{1^{st}-99^{th}}'];
            elseif (variable==2)
                tmp_mean = mean(vr_higher_pct);
                tmp_std = std(vr_higher_pct);
                nt = [mean(dp_higher_pct), std(dp_higher_pct), sum(dp_higher_pct)];
                desp = [wue_name '_{' int2str(pct_best_r_agg(type)) '^{' sup1 '}-99^{th}}'];
            elseif (variable==3)
                tmp_mean = mean(vr_lower_pct);
                tmp_std = std(vr_lower_pct);
                nt = [mean(dp_lower_pct), std(dp_lower_pct), sum(dp_lower_pct)];
                desp = [wue_name '_{1^{st}-' int2str(pct_best_r_agg(type)-1) '^{' sup2 '}}'];
            end
            [tmp, idx] = sort(tmp_mean,'descend');
            col = length(tmp_mean);
            subplot(3, 4, 1+4*(variable-1):3+4*(variable-1))
            hhh = barwitherr(tmp_std(idx(1:col)), tmp_mean(idx(1:col)));
            ylabel({'Predictor importance'}, 'FontSize', 12);
            h = gca;
            h.XTick = 1:col;
            h.YTick = 0:4;
            h.XTickLabel = vr_name(idx(1:col));
            h.TickLabelInterpreter = 'none';

            axis tight
            id = char(erase(id_name(2*variable-1),'"'));
            t = text(0.02,1.22,[id ' ' desp],...
                'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12);
            % most important variable at individual sites
            fff = [];
            if (variable==1)
                tmp = site_vr_agg_all(:,1);
            elseif (variable==2)
                tmp = site_vr_agg_high(:,1);
            elseif (variable==3)
                tmp = site_vr_agg_low(:,1);
            end
            idx_1st = mode(tmp);
            fff = [fff; length(find(tmp==idx_1st))];
            tmp(tmp==idx_1st) = [];
            idx_2nd = mode(tmp);
            fff = [fff; length(find(tmp==idx_2nd))];
            tmp(tmp==idx_2nd) = [];
            idx_3rd = mode(tmp);
            fff = [fff; length(find(tmp==idx_3rd))];
            subplot(3, 4, 4*(variable))
            bar(fff, 0.6)
            ylabel({'Frequency'}, 'FontSize', 12);
            h = gca;
            vr_cal = isfinite(idx_1st)+isfinite(idx_2nd)+isfinite(idx_3rd);
            h.XTick = 1:vr_cal;
            idx_name = [idx_1st, idx_2nd, idx_3rd];
            idx_name(~isfinite(idx_name)) = [];
            h.XTickLabel = vr_name(idx_name);
            h.TickLabelInterpreter = 'none';
            axis tight
            id = char(erase(id_name(2*variable),'"'));
            t = text(0.02,1.22,[id ' ' desp],...
                'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12);
        end
        figname=[path 'RF_VR_WUE_pct_median_ctrl_v2' fig_type '.jpeg'];
        print('-djpeg','-r300',[figname]);
    end
end

%%
% random-forest model variable importance: which factor controls WUE at individual percentiles?
% are the factors contributing to uWUE sensitive to MAH trends?
close all;
clc;
id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", ...
    "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", ...
    "(q)", "(r)", "(s)", "(t)", "(u)", "(v)", "(w)", "(x)"];
for type = 1:1
    if (type==1)
        wue_data = uwue_MAH_pct_ctrl_median_grs;
        path = '../plots_hourly/uWUE_trend_full_v7/uwue_MAH_pct_ctrl_grs/';
        wue_name = '{\it uWUE}';
        sup1 = 'th';
        sup2 = 'th';
        mah_data = wue_hourly_agg;
        pct_threshold_agg = pct_best_r_agg(1);
        col_wue = 3;
    elseif (type==2)
        wue_data = iwue_MAH_pct_ctrl_median_grs;
        path = '../plots_hourly/uWUE_trend_full_v7/iwue_MAH_pct_ctrl_grs/';
        wue_name = '{\it WUE_{ei}}';
        sup1 = 'th';
        sup2 = 'th';
    elseif (type==3)
        wue_data = uwue_MAH_pct_ctrl_median_sum;
        path = '../plots_hourly/uWUE_trend_full_v7/uwue_MAH_pct_ctrl_sum/';
        wue_name = '{\it uWUE}';
        sup1 = 'th';
        sup2 = 'th';
    elseif (type==4)
        wue_data = iwue_MAH_pct_ctrl_median_sum;
        path = '../plots_hourly/uWUE_trend_full_v7/iwue_MAH_pct_ctrl_sum/';
        wue_name = '{\it WUE_{ei}}';
        sup1 = 'rd';
        sup2 = 'nd';
    end
    % calculate MAH trends
    sites = unique(wue_data(:,19));
    site_pos_mah = [];
    site_neg_mah = [];
    for ii = 1:length(sites)
        idx = find(mah_data(:,4)==sites(ii));
        data_tmp = mah_data(idx, :);
        % define growing season as GPP>0 & Tair>0
        idx = find(data_tmp(:, 5)>0 & data_tmp(:, 8)>0);
        grs_data = data_tmp(idx,:);
        yr_range = unique(grs_data(:, 1));
        clear wue_raw MAH_agg
        for yr = 1:length(yr_range) % loop through years
            idx = find(grs_data(:, 1)==yr_range(yr));
            data_site_year = grs_data(idx, :);
            tmp = data_site_year(:, col_wue);
            wue_raw(yr) = median(tmp,'omitnan');
            idx = find(data_site_year(:, col_wue)>=pct_value_agg);
            MAH_agg(yr) = length(idx);
        end
        idx_valid = find(isfinite(wue_raw)==1);
        datain = nan(length(idx_valid),2);
        datain(:,1) = yr_range(idx_valid);
        datain(:,2) = MAH_agg(idx_valid);
        % MAH trends
        [taub tau h sig Z S sigma sen_MAH n senplot CIlower CIupper D Dall C3 nsigma] =...
            ktaub(datain(:, [1, 2]), 0.05);
        if (sen_MAH>0)
            site_pos_mah = [site_pos_mah; sites(ii)];
        elseif (sen_MAH<0)
            site_neg_mah = [site_neg_mah; sites(ii)];
        end
    end
    % group into sites with positive or negative MAH trends
    for group = 1:2
        tmp = wue_data(:,21);
        if (group==1) % positive MAH trends
            sites = site_pos_mah;
            fig_type = '_MAH_p';
        elseif (group==2) % negative MAH trends
            sites = site_neg_mah;
            fig_type = '_MAH_n';
        end        
        r2_lower_pct = [];
        vr_lower_pct = [];
        dp_lower_pct = [];
        r2_higher_pct = [];
        vr_higher_pct = [];
        dp_higher_pct = [];
        r2_all_pct = [];
        vr_all_pct = [];
        dp_all_pct = [];
        site_vr_agg_all = [];
        site_vr_agg_high = [];
        site_vr_agg_low = [];
        for site = 1:length(sites)
            clf;
            siteID = dir(sites(site), 5:10);
            for variable = 1:3
                if (variable==1)
                    idx = find(wue_data(:,19)==sites(site) & wue_data(:,20)<=99);
                    desp = [wue_name '_{1^{st}-99^{th}}'];
                elseif (variable==2)
                    idx = find(wue_data(:,19)==sites(site) & wue_data(:,20)>=pct_best_r_agg(type) & wue_data(:,20)<=99);
                    desp = [wue_name '_{' int2str(pct_best_r_agg(type)) '^{' sup1 '}-99^{th}}'];
                elseif (variable==3)
                    idx = find(wue_data(:,19)==sites(site) & wue_data(:,20)>=1 & wue_data(:,20)<pct_best_r_agg(type));
                    desp = [wue_name '_{1^{st}-' int2str(pct_best_r_agg(type)-1) '^{' sup2 '}}'];
                end
                % build up variable list
                nt = length(idx);
                Tair = wue_data(idx,3);
                VPD = wue_data(idx,5);
                CO2 = wue_data(idx,6);
                SW = wue_data(idx,7);
                USTAR = wue_data(idx,8);
                Wind = wue_data(idx,9);
                PPFD = wue_data(idx,10);
                Tsoil = wue_data(idx,11);
                SWC = wue_data(idx,12);
                G = wue_data(idx,13);
                GPP = wue_data(idx,15);
                ET = wue_data(idx,16);
                WUE = wue_data(idx,1);
                X = table(Tair, VPD, CO2, SW, USTAR, Tsoil, SWC, G, ...
                    WUE);
                t = templateTree('NumPredictorsToSample','all', ...
                    'PredictorSelection', 'curvature', 'Surrogate', 'on');
                rng('default'); % For reproducibility
                Mdl = fitrensemble(X, 'WUE', 'Method', 'bag', 'NumLearningCycles', 100, 'Learners', t);
                imp = oobPermutedPredictorImportance(Mdl);
                [tmp, idx] = sort(imp,'descend');

                % store output for later comparison
                if (variable==3)
                    vr_lower_pct = [vr_lower_pct; imp];
                    r2_lower_pct = [r2_lower_pct; R2];
                    dp_lower_pct = [dp_lower_pct; nt];
                    site_vr_agg_low = [site_vr_agg_low; idx];
                elseif (variable==2)
                    vr_higher_pct = [vr_higher_pct; imp];
                    r2_higher_pct = [r2_higher_pct; R2];
                    dp_higher_pct = [dp_higher_pct; nt];
                    site_vr_agg_high = [site_vr_agg_high; idx];
                elseif (variable==1)
                    vr_all_pct = [vr_all_pct; imp];
                    r2_all_pct = [r2_all_pct; R2];
                    dp_all_pct = [dp_all_pct; nt];
                    site_vr_agg_all = [site_vr_agg_all; idx];
                end
                vr_name = Mdl.PredictorNames;
                vr_name{5} = 'u*';
            end
        end
        % RF variable importance across sites
        clc;
        close all
        id_name = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", ...
            "(i)", "(j)", "(k)", "(l)", "(m)", "(n)", "(o)", "(p)", ...
            "(q)", "(r)", "(s)", "(t)", "(u)", "(v)", "(w)", "(x)"];

        for variable = 1:3
            if (variable==1)
                tmp_mean = mean(vr_all_pct);
                tmp_std = std(vr_all_pct);
                nt = [mean(dp_all_pct), std(dp_all_pct), sum(dp_all_pct)];
                desp = [wue_name '_{1^{st}-99^{th}}'];
            elseif (variable==2)
                tmp_mean = mean(vr_higher_pct);
                tmp_std = std(vr_higher_pct);
                nt = [mean(dp_higher_pct), std(dp_higher_pct), sum(dp_higher_pct)];
                desp = [wue_name '_{' int2str(pct_best_r_agg(type)) '^{' sup1 '}-99^{th}}'];
            elseif (variable==3)
                tmp_mean = mean(vr_lower_pct);
                tmp_std = std(vr_lower_pct);
                nt = [mean(dp_lower_pct), std(dp_lower_pct), sum(dp_lower_pct)];
                desp = [wue_name '_{1^{st}-' int2str(pct_best_r_agg(type)-1) '^{' sup2 '}}'];
            end
            [tmp, idx] = sort(tmp_mean,'descend');
            col = length(tmp_mean);
            subplot(3, 4, 1+4*(variable-1):3+4*(variable-1))
            hhh = barwitherr(tmp_std(idx(1:col)), tmp_mean(idx(1:col)));
            ylabel({'Predictor importance'}, 'FontSize', 12);
            h = gca;
            h.XTick = 1:col;%length(idx);
            h.YTick = 0:4;
            h.XTickLabel = vr_name(idx(1:col));
            h.TickLabelInterpreter = 'none';

            axis tight
            id = char(erase(id_name(2*variable-1),'"'));
            t = text(0.02,1.22,[id ' ' desp],...
                'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12);
            % most important variable at individual sites
            fff = [];
            if (variable==1)
                tmp = site_vr_agg_all(:,1);
            elseif (variable==2)
                tmp = site_vr_agg_high(:,1);
            elseif (variable==3)
                tmp = site_vr_agg_low(:,1);
            end
            idx_1st = mode(tmp);
            fff = [fff; length(find(tmp==idx_1st))];
            tmp(tmp==idx_1st) = [];
            idx_2nd = mode(tmp);
            fff = [fff; length(find(tmp==idx_2nd))];
            tmp(tmp==idx_2nd) = [];
            idx_3rd = mode(tmp);
            fff = [fff; length(find(tmp==idx_3rd))];
            subplot(3, 4, 4*(variable))
            bar(fff, 0.6)
            ylabel({'Frequency'}, 'FontSize', 12);
            h = gca;
            vr_cal = isfinite(idx_1st)+isfinite(idx_2nd)+isfinite(idx_3rd);
            h.XTick = 1:vr_cal;
            idx_name = [idx_1st, idx_2nd, idx_3rd];
            idx_name(~isfinite(idx_name)) = [];
            h.XTickLabel = vr_name(idx_name);
            h.TickLabelInterpreter = 'none';
            axis tight
            id = char(erase(id_name(2*variable),'"'));
            t = text(0.02,1.22,[id ' ' desp],...
                'Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12);
        end
        figname=[path 'RF_VR_WUE_pct_median_ctrl_v2' fig_type '.jpeg'];
        print('-djpeg','-r300',[figname]);
    end
end
