% Code creating array linking month to season for visualisation. Called by
% NPZD_model_forcing.m

% CREATE ARRAYS FOR THE MONTHS (1 to 12) ASSOCIATED WITH EACH DAY IN THE
% SIMULATION
months_label = ["Jan","Feb","Mar","Apr","May","Jun",...
    "Jul","Aug","Sep","Oct","Nov","Dec"];
season_label = ["Summer","Autumn","Winter","Spring"];
if params.hemisphere == "N"                     % If in the Northern Hemisphere
    season_label = circshift(season_label,-2);  % Shift season association to calendar months
end
season_nums = ceil(mod(2:13,12)/3);             % Season index associated with each calendar month
season_nums(season_nums==0) = 4;                % Adjust for season 4
params.month_season = cell(2,12);               % Initialise cell for recording calendar month and season
for indx = 1:12
    params.month_season{1, indx} = months_label(indx);
    params.month_season{2, indx} = season_label(season_nums(indx));
end
cut_off_day = (params.repeat_num-1)*(num_days); % Day at which burn-in ends

params.months_per_day = zeros(1,num_days);          % Initialise array for recording calendar month per day
params.years_per_day = zeros(1, num_days);          % Initialise array for recording calendar year per day
months_days = [31,0,31,30,31,30,31,31,30,31,30,31]; % Number of days in each month
months_nums = circshift(1:12,1-month1);             % Shift to align with first month of data
num_years = ceil(num_days/366);                     % Calculate number of years of data
old_day = 1;
initial_day_of_data = sum(months_days(1:month1-1)) + params.day1; % Calculate the day in the year we are observing
for yr = 0:num_years-1
    if mod(params.year1 + yr,4)==0 && mod(params.year1 + yr,400)==0
        months_days(2)=29;                          % Leap year that is a multiple of 400
    elseif mod(params.year1 + yr,4)==0 && mod(params.year1 + yr,100)~=0
        months_days(2)=29;                          % Every other leap year 
    else
        months_days(2)=28;                          % Not a leap year
    end
    for indx = 1:12
        mnth = months_nums(indx);
        if indx == 1
            new_day = old_day + months_days(mnth)-params.day1;
        else 
            new_day = min(old_day + months_days(mnth)-1, num_days);
        end
        params.months_per_day(old_day:new_day) = mnth;  % Assign appropriate month to day of data
        params.years_per_day(old_day:new_day) = params.year1 + yr; % Assign appropriate year to day of data

        old_day = new_day+1;
        if new_day == num_days
            break
        end
    end
end
params.months_per_day(params.months_per_day==0)=month1; % Correct for 0 entries
params.months_per_day = repmat(params.months_per_day,[1,params.repeat_num]); % repeat for all repeated years
params.years_per_day_burn_in = [];
BIyr = 0;
for BI = 1:params.repeat_num-1
    BIyr = BIyr+1;
    for iii = 0:num_years-1
        BIyr = BIyr+iii;
        if mod(params.year1 + iii,4)==0 && mod(params.year1 + iii,400)==0
            months_days(2)=29;                          % Leap year that is a multiple of 400
        elseif mod(params.year1 + iii,4)==0 && mod(params.year1 + iii,100)~=0
            months_days(2)=29;                          % Every other leap year 
        else
            months_days(2)=28;                          % Not a leap year
        end
        days_in_year = sum(months_days);
        params.years_per_day_burn_in = [params.years_per_day_burn_in, repmat(BIyr,[1,days_in_year])];
    end
end
params.years_per_day = [params.years_per_day_burn_in, params.years_per_day];


% CREATE TICKS FOR PLOTS
months = record.times./(365.25/12);                 % Convert days to fraction of month
season1 = ceil(month1/3);                           % First season index
year_add=0;
months_final = months(months>0);                    % Non-zero months
months_final = months_final(end);                   % Final month to be plotted
tcks = 0:3:ceil(months_final);
num_ticks = length(tcks);                           % Number of seasons to be plotted  
tck_seas = cell(1,num_ticks);
for tck = season1:num_ticks+season1 
    indx = mod(tck,4);                              % Convert to season within 1 year: 1,2,3,4
    indx(indx == 0) = 4;                            % Correct for 0 entries
    if year_add<(params.repeat_num-1)*num_years     % Burn-in period
        if indx == 1
            year_add=year_add+1;
            tck_seas{tck} = "BI year " + num2str(year_add) + " " + ...
                season_label{indx};                     % Tick for the first Month of every year: Year, Season
        else
            tck_seas{tck}=season_label{indx};           % Tick for the first Month of every year: Season
        end
    else
        if indx == 1
            year_add_edited = year_add - (params.repeat_num-1)*num_years;
            tck_seas{tck} = num2str(params.year1 + year_add_edited) + " " + ...
                season_label{indx};                     % Tick for the first Month of every year: Year, Season
            year_add=year_add+1;
        else
            tck_seas{tck}=season_label{indx};           % Tick for the first Month of every year: Season
        end
    end
end 