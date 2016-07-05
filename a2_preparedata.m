%% Clean the workplace
clear all; close all;
a1_init
addpath proc;

%% init work folders 
% TODO: separeate data and results directories ['results' filesep 'filter_data'];
out_dir = ['results' filesep 'filter_data'];

data_dir = 'results';
if ~isdir(out_dir)&&~isempty(out_dir), mkdir(out_dir); end;

%% LOAD RAW: загрузка сырых данных
% First copy xl and do the conversions..
copyfile(fullfile('NKData.xlsx'), fullfile(out_dir,'NKData.xlsx'));
xl2cs(fullfile(data_dir,'NKData.xlsx'),'qdata','quarterly.csv')
xl2cs(fullfile(data_dir,'NKData.xlsx'),'mdata','monthly.csv')


%%  then load the databases:
dq = dbload(fullfile(data_dir,'quarterly.csv'));
dm = dbload(fullfile(data_dir,'monthly.csv'));

%% Prepare quarterly data
% dq.y_gap_f=hpf2(dq.gdp_f);

%% Prepare monthly data
% Inflation
pi_names = {'pi_mom_sa','pi_f_mom_sa'};
%- MoM -> level
dm = dbbatch(dm,'c$0','cumprod(1+dm.$0/100)','stringList=',pi_names);

%% COMBINE: create a quarterly database
db = struct();

% montly data to quarterly data
nameListLast = {'cpi_mom_sa','cpi_f_mom_sa'};
nameListAverage = {'rs_f','rs','s','oil','prem'};

db =dbmerge(...
    dbbatch(dm,'$0','convert(dm.$0,''q'',''method='',''last'')','nameList=',nameListLast, 'fresh',true), ...
    dbbatch(dm,'$0','convert(dm.$0,''q'',''method='', @mean)','nameList=',nameListAverage, 'fresh',true));
% GDP
db.gdp = dq.gdp;
% External output gap
db.gdp_f=dq.gdp_f;

% make logs
nameListLevels = {'cpi_mom_sa','cpi_f_mom_sa','s','oil','gdp','gdp_f'}; % add GDP
nameListRates = {'rs','rs_f','prem'};
%- all except for percentages
db = dbbatch(db,'$0','100*log(db.$0)','nameList=',nameListLevels);
%- log interest rates
db = dbbatch(db,'$0','100*log(1+db.$0/100)','nameList=',nameListRates);

%% SAVE: save the database with processed data
save(fullfile(out_dir,'processed_data.mat'),'db');
% dbsave(fullfile(out_dir,'model_data.csv'),db,inf,'format','%f',...
%   'class',false);
fprintf('finished\n');