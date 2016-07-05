%% Start IRIS if it is not already running:
[~, sysn] = system('hostname');
sysn=strtrim(sysn);

if strcmp(sysn,'WIN-OBBPR9ANCSG')
    try, irisversion;, catch, addpath 'C:\Users\iav\Downloads\IRIS'; irisstartup;, end;
else
    try, irisversion;, catch, addpath 'C:\Users\AIsakov\Downloads\IRIS'; irisstartup;, end;
end
%% DATES: major dates
sfilt = qq(2006,1); % start of filtration
sfcast = qq(2015,4); % start of forecast
efcast = qq(2030,4); % end of forecast
% For reports:
srpt_filt = sfilt; % start of filtration report
erpt_filt = qq(dat2ypf(sfcast)+4,4); % end of filtration report
srpt_fcast = qq(dat2ypf(sfcast)-1,4); % start of forecast report
erpt_fcast = qq(dat2ypf(sfcast)+3,4); % end of forecast report