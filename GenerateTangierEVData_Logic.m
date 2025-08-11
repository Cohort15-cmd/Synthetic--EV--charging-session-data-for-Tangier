function GenerateTangierEVData_Logic()
% Tangier EV charging simulation using a FIXED fleet of 300 vehicles
% (given mix), real stations, operator preference, queuing, and
% plausibility checks. Produces chronologically sorted sessions +
% a daily station report.
%


%% ---------------- Parameters ----------------
outFile   = 'tangier_ev_charging_synthetic.csv';
startDate = datetime(2023,1,1);
endDate   = datetime(2024,12,31);
cityName  = 'Tangier';

% Time-of-day shape for plug-in times (not volume)
lambda0_time_shape = 2.5;    % larger -> sharper AM/PM peaks

% Queuing (max willing to wait)
maxQueueAC = 15;   % minutes
maxQueueDC = 20;   % minutes

% Daily driving intensity (km): scale factors
WD_scale = 1.0;      % weekdays
WE_scale = 0.8;      % weekends (less driving on average)
SUMMER_scale = 1.1;  % Jun–Aug a bit higher

%% ---------------- Stations (verified list) ----------------
% {Name, Lat, Lon, Max_kW, ConnectorTypes, StallCount, Operator, SupportsNonTesla, ConnectorCounts}
stations = {
    'Tesla_Supercharger_Hilton_AlHouara', 35.666741, -5.965834, 150, {'Tesla','CCS2'}, 4, 'Tesla', false, struct('Tesla',4,'CCS2',4,'CHAdeMO',0,'Type2',0);
    'FastVolt_Afriquia_Mnar',            35.8020831, -5.7436656, 100, {'CCS2','CHAdeMO','Type2'}, 3, 'FastVolt', true, struct('Tesla',0,'CCS2',1,'CHAdeMO',1,'Type2',2);
    'FastVolt_Afriquia_TFZ',             35.726800, -5.882000, 50,  {'CCS2','CHAdeMO','Type2'}, 3, 'FastVolt', true, struct('Tesla',0,'CCS2',1,'CHAdeMO',1,'Type2',1);
    'FastVolt_Afriquia_Golf',            35.7679901, -5.8596349, 22, {'Type2'}, 2, 'FastVolt', true, struct('Tesla',0,'CCS2',0,'CHAdeMO',0,'Type2',2);
    'FastVolt_Marjane_Tangier',          35.780100, -5.813600, 22, {'Type2'}, 3, 'FastVolt', true, struct('Tesla',0,'CCS2',0,'CHAdeMO',0,'Type2',3);
    'TotalEnergies_Relais_Tanger',       35.738923, -5.852115, 22, {'Type2'}, 1, 'TotalEnergies', true, struct('Tesla',0,'CCS2',0,'CHAdeMO',0,'Type2',1);
    'Tanger_City_Mall',                  35.780650, -5.799300, 22, {'Type2'}, 2, 'Public', true, struct('Tesla',0,'CCS2',0,'CHAdeMO',0,'Type2',2);
    'Kenzi_Solazur_Hotel',               35.779500, -5.790800, 11, {'Type2'}, 2, 'Hotel', true, struct('Tesla',0,'CCS2',0,'CHAdeMO',0,'Type2',2);
    'Royal_Golf_Tanger',                 35.744200, -5.849300, 11, {'Type2'}, 2, 'Recreation', true, struct('Tesla',0,'CCS2',0,'CHAdeMO',0,'Type2',2);
    'Royal_Tulip_CityCenter',            35.776500, -5.834000, 11, {'Type2'}, 1, 'Hotel', true, struct('Tesla',0,'CCS2',0,'CHAdeMO',0,'Type2',1);
    'Barcelo_Occidental_Tanger',         35.779000, -5.817000, 11, {'Type2'}, 1, 'Hotel', true, struct('Tesla',0,'CCS2',0,'CHAdeMO',0,'Type2',1)
};

%% ---------------- Fleet (exact 300 units) ----------------
% Build the 300-vehicle fleet with per-vehicle specs
fleet = build_fleet();     % table: vehicle_id, brand, model, battKwh, conn, PvehAC, PvehDC
N = height(fleet);

% Initial SOC state for each vehicle (reflects mixed prior usage)
socNow = 60 + 30*rand(N,1);  % 60–90%

%% ---------------- Runtime station state (queuing) ----------------
numStations = size(stations,1);
stallCountsVec = zeros(numStations,1);
for si=1:numStations
    sc = stations{si,6}; if isnan(sc) || sc<=0, sc = 1; end
    stallCountsVec(si) = sc;
end
stationNextFree = cell(numStations,1);
for si=1:numStations
    stationNextFree{si} = repmat(datetime(-Inf,'ConvertFrom','posixtime'), stallCountsVec(si), 1);
end

%% ---------------- Simulation loop ----------------
allDays = startDate:endDate;
rows   = {};

% Persistent anomaly state for climate model
prevAnom = NaN;

for d = 1:numel(allDays)
    baseDay = allDays(d);

    % --- Daily temperature state (uses persistence across days)
    m = month(baseDay);
    [TminDay, TmaxDay, prevAnom] = daily_temp_Tangier(m, prevAnom);

    % Driving scale
    wdName   = char(day(baseDay,'longname'));
    wdScale  = (ismember(wdName, {'Saturday','Sunday'})) * WE_scale + ...
               (~ismember(wdName, {'Saturday','Sunday'})) * WD_scale;
    if month(baseDay)>=6 && month(baseDay)<=8, wdScale = wdScale * SUMMER_scale; end

    % --- For each vehicle: consume energy by daily driving, maybe schedule charge
    todaysRequests = [];
    for vi = 1:N
        B = fleet.battKwh(vi);
        D = daily_distance_km() * wdScale;
        % approximate energy based on daily mean temperature (midpoint)
        Tmid = 0.5*(TminDay + TmaxDay);
        Etrip = D * energy_per_km(Tmid);
        dSOC  = 100*Etrip/B;
        socNow(vi) = max(0, socNow(vi) - dSOC + randn*0.5);

        % Decide if this vehicle will attempt to charge today based on range
        range_km = (B * socNow(vi)/100) / energy_per_km(Tmid);
        Rth_AC = 40;  k_AC = 20;
        Rth_DC = 120; k_DC = 25;
        p_ACneed = 1 ./ (1 + exp((range_km - Rth_AC)/k_AC));
        p_DCneed = 1 ./ (1 + exp((range_km - Rth_DC)/k_DC));
        p_need   = max(p_ACneed, 0.25*p_DCneed);

        if rand < p_need
            ts_req = sample_one_time(baseDay, lambda0_time_shape);
            todaysRequests(end+1,:) = [vi, posixtime(ts_req)]; 
        end
    end

    % Sort today's requests by time and process
    if ~isempty(todaysRequests)
        todaysRequests = sortrows(todaysRequests, 2);
        for r = 1:size(todaysRequests,1)
            vi    = todaysRequests(r,1);
            tsArr = datetime(todaysRequests(r,2), 'ConvertFrom','posixtime');

            if socNow(vi) > 95, continue; end

            evName  = fleet.model{vi};
            vehConn = fleet.conn{vi};
            B       = fleet.battKwh(vi);
            PvehDC_cap = fleet.PvehDC(vi);

            stIdx = choose_station(evName, vehConn, stations);
            stationID = stations{stIdx,1};
            lat = stations{stIdx,2};
            lon = stations{stIdx,3};
            siteMaxKW = stations{stIdx,4};
            connList  = stations{stIdx,5};
            numStalls = stations{stIdx,6}; if isnan(numStalls) || numStalls<=0, numStalls = 1; end
            supportsNonTS = stations{stIdx,8};
            connCounts = stations{stIdx,9};

            hasDC = any(ismember(connList, {'CCS2','CHAdeMO','Tesla'}));
            hasAC = any(ismember(connList, {'Type2'}));

            canUseDC = false;
            if PvehDC_cap > 0
                switch vehConn
                    case 'CCS2'
                        canUseDC = (connCounts.CCS2 > 0) || ((connCounts.Tesla > 0) && supportsNonTS);
                    case 'CHAdeMO'
                        canUseDC = (connCounts.CHAdeMO > 0);
                    case 'Tesla'
                        canUseDC = (connCounts.Tesla > 0) || (connCounts.CCS2 > 0);
                    otherwise
                        canUseDC = false;
                end
            end

            % Decide DC vs AC
            T_arr = temp_at_time(tsArr, TminDay, TmaxDay);
            range_km = (B * socNow(vi)/100) / energy_per_km(T_arr);
            Rth_AC = 40;  k_AC = 20;
            Rth_DC = 120; k_DC = 25;
            p_AC = 1 ./ (1 + exp((range_km - Rth_AC)/k_AC));
            p_DC = 1 ./ (1 + exp((range_km - Rth_DC)/k_DC));

            if hasDC && canUseDC
                preferDC = (range_km < 150) || (socNow(vi) < 60);
                if preferDC || ~hasAC
                    isDC = rand < max(p_DC,0.5);
                else
                    isDC = rand < 0.3;
                end
            else
                isDC = false;
            end
            if ~isDC && ~hasAC && hasDC && canUseDC
                isDC = true;
            end

            p_charge = isDC*p_DC + (~isDC)*p_AC;
            if rand > p_charge
                continue;
            end

            siteMaxKW_orig = siteMaxKW;
            if isDC
                if siteMaxKW >= 130
                    chargerType = 'DC_150kW'; tierKW = 150;
                elseif siteMaxKW >= 90
                    chargerType = 'DC_100kW'; tierKW = 100;
                elseif siteMaxKW >= 50
                    chargerType = 'DC_50kW';  tierKW = 50;
                else
                    chargerType = sprintf('DC_%dkW', siteMaxKW); tierKW = siteMaxKW;
                end
                minDelta = 10; targetSet = [70 80 85 90]; targetProb = [0.1 0.5 0.3 0.1];
            else
                if siteMaxKW >= 22
                    chargerType = 'AC_22kW';  tierKW = 22;
                elseif siteMaxKW >= 11
                    chargerType = 'AC_11kW';  tierKW = 11;
                else
                    chargerType = sprintf('AC_%dkW', siteMaxKW); tierKW = siteMaxKW;
                end
                minDelta = 5;  targetSet = [80 90 100]; targetProb = [0.25 0.5 0.25];
            end

            socArr = socNow(vi);
            socDep = weighted_pick(targetSet, targetProb);
            socDep = max(socDep, socArr + minDelta);
            socDep = min(socDep, 100);

            [durMin, energyKwh] = duration_from_soc(B, socArr, socDep, T_arr, isDC, evName, tierKW);

            avgP = (durMin>0) * (energyKwh / (durMin/60));
            if ~plausible_session(isDC, socArr, socDep, durMin, avgP, tierKW)
                continue;
            end

            [nextWaitMin, stallIdx, okQueue] = queue_decision(stationNextFree{stIdx}, tsArr, isDC, durMin, maxQueueAC, maxQueueDC);
            if ~okQueue
                continue;
            end
            plugTs = tsArr + minutes(nextWaitMin);
            depTs  = plugTs + minutes(durMin);
            stationNextFree{stIdx}(stallIdx) = depTs;

            socNow(vi) = socDep - 0.5 + randn*0.25;
            socNow(vi) = min(max(socNow(vi), 0), 100);

            dayName = char(day(tsArr,'longname'));
            isWeekend = ismember(dayName, {'Saturday','Sunday'});

            rows(end+1,:) = {
                sprintf('S%05d', size(rows,1)+1), ...
                fleet.vehicle_id{vi}, ...
                stationID, cityName, lat, lon, tsArr, dayName, logical(isWeekend), ...
                chargerType, evName, B, round(socArr,1), round(socDep,1), round(T_arr,1), ...
                round(durMin,1), round(energyKwh,2), siteMaxKW_orig, tierKW, ...
                numStalls, stations{stIdx,7}, stations{stIdx,8}, jsonencode(connCounts), strjoin(connList,'/'), ...
                round(nextWaitMin,1), plugTs, depTs, stallIdx}; %#ok<AGROW>
        end
    end
end

%% ---------------- Write CSV ----------------
if isempty(rows)
    warning('No sessions generated.');
    return;
end
T = cell2table(rows, 'VariableNames', ...
    {'session_id','vehicle_id','station_id','city','lat','lon','timestamp_arrival', ...
     'day_of_week','is_weekend','charger_type','ev_model','battery_kwh', ...
     'arrival_soc_pct','departure_soc_pct','ambient_temp_c', ...
     'duration_min','energy_kwh','station_max_kw','site_power_tier_kw','num_stalls','operator','supports_non_tesla','port_counts','station_connectors', ...
     'queue_wait_min','timestamp_plug_in','timestamp_departure','stall_index'});
T = sortrows(T, 'timestamp_arrival');

writetable(T, outFile);
fprintf('Saved %d sessions to %s', height(T), outFile);

%% ===================== Daily Station Report =====================
reportRows = {};
if height(T) > 0
    stallMap = containers.Map();
    for si = 1:size(stations,1)
        stallMap(stations{si,1}) = stations{si,6};
    end
    T.day = dateshift(T.timestamp_arrival, 'start', 'day');
    [grp, Gstation, Gday] = findgroups(T.station_id, T.day);
    sessCount = splitapply(@numel, T.session_id, grp);
    energyKWh = splitapply(@sum, T.energy_kwh, grp);
    stallMins = splitapply(@sum, T.duration_min, grp);
    avgWait   = splitapply(@mean, T.queue_wait_min, grp);
    maxWait   = splitapply(@max,  T.queue_wait_min, grp);
    peakConc = zeros(size(sessCount));
    for gi = 1:numel(sessCount)
        st = Gstation{gi}; d0 = Gday(gi);
        mask = strcmp(T.station_id, st) & (T.day == d0);
        starts = T.timestamp_plug_in(mask);
        ends   = T.timestamp_departure(mask);
        n = sum(mask);
        if n==0, peakConc(gi)=0; continue; end
        deltas  = [ones(n,1); -ones(n,1)];
        [~, idx] = sort([starts; ends]);
        deltas = deltas(idx);
        c = 0; m = 0;
        for k=1:numel(deltas)
            c = c + deltas(k);
            if c>m, m=c; end
        end
        peakConc(gi) = m;
    end
    for gi = 1:numel(sessCount)
        st = Gstation{gi}; d0 = Gday(gi);
        stalls = 1; if isKey(stallMap, st), stalls = stallMap(st); end
        util = 0; if stalls>0, util = (stallMins(gi) / (stalls*1440)) * 100; end
        busy = util >= 60;
        reportRows(end+1,:) = {st, d0, stalls, sessCount(gi), energyKWh(gi), stallMins(gi), util, peakConc(gi), avgWait(gi), maxWait(gi), busy};
    end
    R = cell2table(reportRows, 'VariableNames', ...
        {'station_id','day','num_stalls','sessions','energy_kwh','stall_minutes','utilization_pct','peak_concurrent','avg_queue_min','max_queue_min','busy_flag'});
    writetable(R, 'station_daily_report.csv');
    fprintf('Saved station daily report: %d rows -> station_daily_report.csv', height(R));
else
    fprintf('No sessions -> station_daily_report.csv not created.');
end

end

%% ===================== Fleet builder & specs =====================
function F = build_fleet()
% Construct the 300-EV fleet from the provided distribution.
% Returns table with columns:
% vehicle_id, brand, model, battKwh, conn, PvehAC, PvehDC

specs = { ... % brand, model, units
'Dacia','Spring',        122; ...
'BYD','Dolphin',          23; ...
'BYD','Atto 3',           10; ...
'Volvo','XC40 Recharge',  23; ...
'Volvo','C40 Recharge',    6; ...
'DFSK','Seres 3',         13; ...
'DFSK','EC35 Van',         8; ...
'Citroën','ë-C4',         15; ...
'Citroën','ë-Berlingo',    6; ...
'Mercedes','EQB',          6; ...
'Mercedes','EQA',          6; ...
'Fiat','500e',            11; ...
'Renault','Megane E-Tech', 6; ...
'Renault','Kangoo E-Tech', 3; ...
'Hyundai','Kona Electric', 6; ...
'Hyundai','Ioniq 5',       2; ...
'Audi','Q4 e-tron',        8; ...
'Kia','EV6',               5; ...
'Kia','Niro EV',           2; ...
'Opel','Mokka-e',          4; ...
'Opel','Corsa-e',          2; ...
'MG','MG4',                4; ...
'MG','MG ZS EV',           2; ...
'Geely','Geometry C',      4; ...
'BMW','iX3',               2; ...
'BMW','i4',                2; ...
'Tesla','Model 3',         5; ...
'Tesla','Model Y',         4  ...
};

rows = {};
vid = 1;
for i=1:size(specs,1)
    brand = specs{i,1}; model = specs{i,2}; units = specs{i,3};
    for k=1:units
        name = sprintf('%s_%s', sanitize(brand), sanitize(model));
        [B, conn, Pac, Pdc] = vehicle_specs(name);
        rows(end+1,:) = {sprintf('V%03d', vid), brand, name, B, conn, Pac, Pdc};
        vid = vid + 1;
    end
end
F = cell2table(rows, 'VariableNames', {'vehicle_id','brand','model','battKwh','conn','PvehAC','PvehDC'});
end

function s = sanitize(x)
s = regexprep(x, '[^A-Za-z0-9]+', '_');
end

function [B, conn, Pac, Pdc] = vehicle_specs(model)
switch model
    case 'Dacia_Spring',        B=27;   conn='CCS2'; Pac=7.4;  Pdc=30;
    case 'BYD_Dolphin',         B=60;   conn='CCS2'; Pac=11;   Pdc=90;
    case 'BYD_Atto_3',          B=60;   conn='CCS2'; Pac=11;   Pdc=88;
    case 'Volvo_XC40_Recharge', B=75;   conn='CCS2'; Pac=11;   Pdc=150;
    case 'Volvo_C40_Recharge',  B=75;   conn='CCS2'; Pac=11;   Pdc=150;
    case 'DFSK_Seres_3',        B=52;   conn='CCS2'; Pac=6.6;  Pdc=90;
    case 'DFSK_EC35_Van',       B=38;   conn='CCS2'; Pac=6.6;  Pdc=40;
    case 'Citro_n___C4',        B=50;   conn='CCS2'; Pac=11;   Pdc=100;
    case 'Citro_n___Berlingo',  B=50;   conn='CCS2'; Pac=11;   Pdc=100;
    case 'Mercedes_EQB',        B=66;   conn='CCS2'; Pac=11;   Pdc=120;
    case 'Mercedes_EQA',        B=66;   conn='CCS2'; Pac=11;   Pdc=120;
    case 'Fiat_500e',           B=42;   conn='CCS2'; Pac=11;   Pdc=85;
    case 'Renault_Megane_E_Tech', B=60; conn='CCS2'; Pac=22;   Pdc=130;
    case 'Renault_Kangoo_E_Tech', B=45; conn='CCS2'; Pac=22;   Pdc=80;
    case 'Hyundai_Kona_Electric', B=64; conn='CCS2'; Pac=11;   Pdc=77;
    case 'Hyundai_Ioniq_5',     B=73;   conn='CCS2'; Pac=11;   Pdc=220;
    case 'Audi_Q4_e_tron',      B=77;   conn='CCS2'; Pac=11;   Pdc=125;
    case 'Kia_EV6',             B=77;   conn='CCS2'; Pac=11;   Pdc=220;
    case 'Kia_Niro_EV',         B=65;   conn='CCS2'; Pac=11;   Pdc=77;
    case 'Opel_Mokka_e',        B=50;   conn='CCS2'; Pac=11;   Pdc=100;
    case 'Opel_Corsa_e',        B=50;   conn='CCS2'; Pac=11;   Pdc=100;
    case 'MG_MG4',              B=64;   conn='CCS2'; Pac=11;   Pdc=140;
    case 'MG_MG_ZS_EV',         B=44.5; conn='CCS2'; Pac=11;  Pdc=75;
    case 'Geely_Geometry_C',    B=70;   conn='CCS2'; Pac=11;   Pdc=75;
    case 'BMW_iX3',             B=74;   conn='CCS2'; Pac=11;   Pdc=150;
    case 'BMW_i4',              B=80;   conn='CCS2'; Pac=11;   Pdc=200;
    case 'Tesla_Model_3',       B=60;   conn='Tesla'; Pac=11;  Pdc=170;
    case 'Tesla_Model_Y',       B=75;   conn='Tesla'; Pac=11;  Pdc=200;
    otherwise,                  B=55;   conn='CCS2'; Pac=11;   Pdc=100;
end
end

%% ===================== Time & demand helpers =====================
function ts = sample_one_time(dayDate, lambda0_time_shape)
    h = (0:23)';
    gauss = @(x,mu,s) exp(-0.5*((x-mu)/s).^2);
    intraday = 0.3*gauss(h,9,1.2) + 0.5*gauss(h,19.5,1.8) + 0.2/24;
    if nargin < 2 || isempty(lambda0_time_shape), lambda0_time_shape = 1; end
    intraday = intraday .^ lambda0_time_shape;
    intraday = intraday / sum(intraday);
    c = cumsum(intraday); r = rand; hh = find(r<=c,1,'first')-1;
    mm = rand*60; ts = dayDate + hours(hh) + minutes(mm);
end

function D = daily_distance_km()
r = rand;
if r < 0.45
    mu = log(12); sigma = 0.5;
elseif r < 0.85
    mu = log(35); sigma = 0.5;
else
    mu = log(90); sigma = 0.6;
end
D = lognrnd(mu, sigma);
D = min(D, 250);
end

%% ===================== Physics & queueing =====================
function eta = eff_from_T(T)
    eta = 0.94 - 0.0008*(25 - T).^2;
    eta = min(max(eta, 0.85), 0.98);
end

function ekm = energy_per_km(T)
    e0 = 0.16; a = 0.008; b = 0.004;
    ekm = e0 * (1 + a*max(0,20-T) + b*max(0,T-30));
end

function dT = derate_T(T)
    dT = 1 - 0.01*max(0, T - 35);
    dT = max(dT, 0.85);
end

function tau = taper_soc(SOC, isDC)
    if SOC <= 80
        tau = 1.0;
    else
        beta = isDC*0.7 + (~isDC)*0.5;
        tau = 1 - beta * (SOC - 80)/20;
        tau = max(tau, 0.2);
    end
end

function P = power_at_SOC(SOC, T, isDC, evName, siteMaxKW)
    Psite = siteMaxKW;
    [Pac,Pdc] = vehicle_caps_quick(evName);
    Pveh  = isDC*Pdc + (~isDC)*Pac;
    Pmax  = min(Psite, Pveh);
    if Pmax <= 0, Pmax = 6; end
    P = Pmax * derate_T(T) * taper_soc(SOC, isDC);
    P = max(P, 3);
end

function [durMin, energyKwh] = duration_from_soc(B, SOCa, SOCt, T, isDC, evName, siteMaxKW)
    SOCa = max(0, min(100, SOCa));
    SOCt = max(SOCa, min(100, SOCt));
    eta = eff_from_T(T);
    socs = ceil(SOCa):1:floor(SOCt);
    if isempty(socs), durMin = 0; energyKwh = 0; return; end
    dtMin = 0; E = 0;
    for s = socs
        P = power_at_SOC(s, T, isDC, evName, siteMaxKW);
        dE = (B/100)/eta;
        dT = 60 * dE / P;
        dtMin = dtMin + dT; E = E + dE;
    end
    durMin = min(max(dtMin, 5), 300);
    energyKwh = E;
end

function [Pac,Pdc] = vehicle_caps_quick(evName)
switch evName
    case 'Dacia_Spring',        Pac=7.4;  Pdc=30;
    case 'BYD_Dolphin',         Pac=11;   Pdc=90;
    case 'BYD_Atto_3',          Pac=11;   Pdc=88;
    case 'Volvo_XC40_Recharge', Pac=11;   Pdc=150;
    case 'Volvo_C40_Recharge',  Pac=11;   Pdc=150;
    case 'DFSK_Seres_3',        Pac=6.6;  Pdc=90;
    case 'DFSK_EC35_Van',       Pac=6.6;  Pdc=40;
    case 'Citro_n___C4',        Pac=11;   Pdc=100;
    case 'Citro_n___Berlingo',  Pac=11;   Pdc=100;
    case 'Mercedes_EQB',        Pac=11;   Pdc=120;
    case 'Mercedes_EQA',        Pac=11;   Pdc=120;
    case 'Fiat_500e',           Pac=11;   Pdc=85;
    case 'Renault_Megane_E_Tech', Pac=22; Pdc=130;
    case 'Renault_Kangoo_E_Tech', Pac=22; Pdc=80;
    case 'Hyundai_Kona_Electric', Pac=11; Pdc=77;
    case 'Hyundai_Ioniq_5',     Pac=11;   Pdc=220;
    case 'Audi_Q4_e_tron',      Pac=11;   Pdc=125;
    case 'Kia_EV6',             Pac=11;   Pdc=220;
    case 'Kia_Niro_EV',         Pac=11;   Pdc=77;
    case 'Opel_Mokka_e',        Pac=11;   Pdc=100;
    case 'Opel_Corsa_e',        Pac=11;   Pdc=100;
    case 'MG_MG4',              Pac=11;   Pdc=140;
    case 'MG_MG_ZS_EV',         Pac=11;   Pdc=75;
    case 'Geely_Geometry_C',    Pac=11;   Pdc=75;
    case 'BMW_iX3',             Pac=11;   Pdc=150;
    case 'BMW_i4',              Pac=11;   Pdc=200;
    case 'Tesla_Model_3',       Pac=11;   Pdc=170;
    case 'Tesla_Model_Y',       Pac=11;   Pdc=200;
    otherwise,                  Pac=11;   Pdc=100;
end
end

function ok = plausible_session(isDC, socArr, socDep, durMin, avgP, tierKW)
    if socArr >= 98 && (socDep - socArr) < 3, ok=false; return; end
    if isDC && socArr > 85 && (socDep - socArr) < 8 && durMin > 25, ok=false; return; end
    if ~isfinite(avgP) || avgP <= 1.5 || avgP > tierKW*1.15, ok=false; return; end
    if ~isDC && durMin > 480, ok=false; return; end
    ok = true;
end

function [waitMin, stallIdx, ok] = queue_decision(nextFreeVec, tsArr, isDC, durMin, maxQac, maxQdc)
    availIdx = find(nextFreeVec <= tsArr, 1, 'first');
    if ~isempty(availIdx)
        waitMin = 0; stallIdx = availIdx; ok = true; return;
    end
    [soonest, idx] = min(nextFreeVec);
    waitMin = minutes(soonest - tsArr);
    base = isDC * maxQdc + (~isDC) * maxQac;
    elastic = 0.8 + 0.002 * min(durMin, 180);
    thr = base * elastic;
    ok  = waitMin <= thr;
    stallIdx = idx;
end

function stIdx = choose_station(evName, vehConn, stations)
    n = size(stations,1); w = zeros(n,1) + 1e-6;
    isTesla = startsWith(evName,'Tesla_');
    for i=1:n
        siteMax = stations{i,4}; conns = stations{i,5}; op = stations{i,7}; allowNonT = stations{i,8}; counts = stations{i,9};
        hasType2 = any(strcmp(conns,'Type2')) && counts.Type2>0;
        hasCCS2  = any(strcmp(conns,'CCS2'))  && counts.CCS2>0;
        hasCHA   = any(strcmp(conns,'CHAdeMO'))&& counts.CHAdeMO>0;
        hasTesla = any(strcmp(conns,'Tesla'))  && counts.Tesla>0;
        compatAny = hasType2 || ...
            (strcmp(vehConn,'CCS2')    && (hasCCS2 || (hasTesla && allowNonT))) || ...
            (strcmp(vehConn,'CHAdeMO') && hasCHA) || ...
            (strcmp(vehConn,'Tesla')   && (hasTesla || hasCCS2));
        if ~compatAny, w(i)=0; continue; end
        wi = 1.0;
        wi = wi * min(2.0, max(0.6, siteMax/22));
        if isTesla
            if strcmp(op,'Tesla'), wi = wi * 4.0; else, wi = wi * 0.6; end
        else
            if strcmp(op,'Tesla'), wi = wi * 0.4; end
            if strcmp(op,'FastVolt'), wi = wi * 1.5; end
        end
        dcBonus = 1.0;
        if strcmp(vehConn,'CCS2') && (hasCCS2 || (hasTesla && allowNonT)), dcBonus = 1.3; end
        if strcmp(vehConn,'CHAdeMO') && hasCHA, dcBonus = 1.3; end
        if strcmp(vehConn,'Tesla') && (hasTesla || hasCCS2), dcBonus = 1.3; end
        w(i) = wi * dcBonus;
    end
    if all(w<=0), stIdx = randi(n); return; end
    c = cumsum(w/sum(w)); r = rand; stIdx = find(r<=c,1,'first');
end

%% ===================== Utility picker =====================
function v = weighted_pick(values, probs)
    probs = probs(:) / sum(probs);
    c = cumsum(probs); r = rand;
    idx = find(r <= c, 1, 'first');
    v = values(idx);
end

%% ===================== Climate helpers (Tangier) =====================
function [tminClim, tmaxClim] = tangier_monthly_normals()
    tmaxClim = [16 16 17 19 21 24 27 28 26 23 19 17];
    tminClim = [ 9 10 10 11 14 16 19 19 18 15 12 10];
end

function [Tmin, Tmax, anomOut] = daily_temp_Tangier(monthIdx, anomIn)
    [tminClim, tmaxClim] = tangier_monthly_normals();
    phi = 0.80;
    sig = [2.2 2.0 2.0 1.8 1.6 1.4 1.3 1.3 1.5 1.8 2.0 2.2];
    if isempty(anomIn) || ~isfinite(anomIn)
        anom = 0;
    else
        anom = phi*anomIn;
    end
    anom = anom + randn*sig(monthIdx);
    if ismember(monthIdx, 6:9) && rand < 0.05
        anom = anom + (2 + 2*rand);
    end
    Tmin = tminClim(monthIdx) + anom + 0.5*randn;
    Tmax = tmaxClim(monthIdx) + anom + 0.5*randn;
    if Tmin > Tmax - 2
        Tmin = Tmax - 2;
    end
    anomOut = anom;
end

function T = temp_at_time(ts, Tmin, Tmax)
    h = hour(ts) + minute(ts)/60;
    T = Tmin + 0.5*(Tmax - Tmin)*(1 - cos(2*pi*(h - 5)/24));
end
