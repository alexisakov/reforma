% Процедура прогнозирования
 
% !!! Блочная структура файла исключительно для удобства навигации по нему,
% запускать его нужно только целиком!!!
% Для обсчета отдельных сценариев надо использовать соответствующие 
% флажки (sim_bl, sim_sc1 и т.п.)

%% очистка рабочего пространства
clear all; close all;
%% выбор симулируемых блоков
% Название сценария
hdstr = 'нефть - 50 в среднем в 2015 г., вклад фондов = 0,5 стандартного предположения,новая инфляция 11,4 в 2014г., фунд. курс 60, без коррекции инфляции';

% базовый сценарий
sim_bl = true;%false

% альт.сценарий 1
sim_sc1 = false;%true
% альт.сценарий 2
sim_sc2 = false;%true;%

% анализ изменений
nlys = false; %true;%  %пока не используется

%% рабочие каталоги
addpath proc;
in_dir = 'res_filter_PIEFnew\';
% out_dir = 'res_fcast_60_pubsec_60_3Y_nostim\';
out_dir = 'Baseline(oil60)\';

if ~isdir(out_dir)&&~isempty(out_dir), mkdir(out_dir); end;

copyfile([mfilename '.m'], [out_dir mfilename '.m']);
% 
load([in_dir 'bl_fcast_ntf_oil40_bsk80.mat'],'fdb');
forecastBASELINE=fdb;
clear fdb;

load([in_dir 'bl_fcast_ntf_tr.mat'],'fdb');
fdb2=fdb;
clear fdb;


%% переключатели предпочтений
fl_ntf = true; % false; % используем ли краткосрочные прогнозы
if fl_ntf, ntfsfx = '_ntf'; else ntfsfx = ''; end;

%% загрузка модели
% m = model('qpm_premSS08_resRSn.mod','linear',true);
m = model('qpm_pieexp_rsn_crisis.mod','linear',true);

load([in_dir 'model_filtration.mat'],'mfilt');

%% изменение базовой параметризации
% стац. состояние
m.PIE_TAR_SS = 100*log(4*0.01+1);% для цели
m.DOT_LTOT_SS = 0; % для условий торговли

% параметры - техническое
m.psi = 0; %AR для PIE_TAR 
m.c1_pie_utils = 0; %AR для PIE_UTL
m.phi4 = .0; %AR для RRF_GAP
% m.chi    = 0.95; %AR для PREM
m.std_RES_OBS_EMBI = 0;% чтобы уж наверняка)

m.w_pie_mp = 0; %просто чтоб не забыть, что оно туттаг

m.w_pie = .9;

m.const_y_growth = 2.0;

%параметры - сценарное
% m.beta4=0;%.03; %волатильный курс не влияет на выпуск
m.gamma2 = 1.9; %бодрее реакция на инфляцию (сценарное мартовского уточнения)
m.gamma3 = 0; % нет реакции на выпуск (сценарное мартовского уточнения)

%% решение модели
m = solve(m);
m = sstate(m);

%% проверка даты фильтрации
res1 = dir([in_dir 'filtered_dbase.mat']);
res2 = dir('a2_filtration.m');
res3 = dir(get(m,'file'));
pause on;
if isempty(res1)
    error(['Без фильтрации не прогнозируем! Нету результатов фильтрации',...
    ' в указанной папке...']);
    pause(2); return
end;
fprintf(['\n\tПроверь себя!\n\tВремя последней фильтрации: %s.',...
    '\n\tПоследние изменения файла [filtration.m]: %s.',...
    '\n\tПоследнее изменение модельного файла: %s.\n\n'],...
    res1.date,res2.date,res3.date);
pause(1);

%% начальные условия и загрузка данных

load('res_dates/dates.mat');

load([in_dir 'filtered_dbase.mat'],'filtdb');
load([in_dir 'obs_data.mat'],'obsdb');
initdb.mean = dbclip(filtdb.mean,sfilt:sfcast-1);
initdb.std = dbclip(filtdb.std,sfilt:sfcast-1);
initdb.mse = filtdb.mse{sfilt:sfcast-1,:,:};

% tmpdb = dbload('hntf_2014q4.csv'); %!!!чтобы учесть изменения в данных, не трогая фильтрацию

%% краткосрочные уточняющие прогнозы

cond = initdb.mean; % база данных для прогноза
obsdb = loadstruct([in_dir 'obs_data.mat']); % база данных, из которой 
obsdb = obsdb.obsdb;                         % будем подцеплять КСП

% % % % % % % % tmpdb = dbload('hntf_2014q3.csv');

p = plan(m, sfilt:efcast);

% КСП по компонентам инфляции

cond.OBS_PIE_P = obsdb.OBS_PIE_P; %продовольствие
cond.OBS_PIE_P(sfcast:sfcast+2) = 100*log(1+0.01*([28.8472 46.9025  29.0173])); %[25.5 21.7 19.4]  %25.584 34.147 24.581  28.16 40.33 -> 26.13 40.08 -> 25.584 81.713 56.883
p = exogenize(p, 'OBS_PIE_P', sfcast:sfcast+2);
p = endogenize(p, 'RES_PIE_P', sfcast:sfcast+2);
% cond.RES_PIE_P(sfcast+3:sfcast+3)= [-(0.5*8.0661)]; %ликвидация эффекта от курсовой динамики (накопившейся в 2014Q4-2015Q1
% cond.RES_PIE_P(sfcast+4:sfcast+6) = arf(tseries(sfcast+3,-(0.5*8.0661)),[1 -.8],0,sfcast+4:sfcast+6);

cond.OBS_PIE_NP = obsdb.OBS_PIE_NP; %непрод.тов.
cond.OBS_PIE_NP(sfcast:sfcast+2) = 100*log(1+0.01*([14.7496 23.1089 15.5921])); %14.839 20.239 15.592  0.97 24.53 -> 14.82 24.57 -> 14.839 15.543 29.56
p = exogenize(p, 'OBS_PIE_NP', sfcast:sfcast+2);
p = endogenize(p, 'RES_PIE_NP', sfcast:sfcast+2);
% cond.RES_PIE_NP(sfcast+3:sfcast+3)= 0.5*[-10.7547]; %ликвидация эффекта от курсовой динамики (накопившейся в 2014Q4-2015Q1
% cond.RES_PIE_NP(sfcast+4:sfcast+6) = arf(tseries(sfcast+3,cond.RES_PIE_NP(sfcast+3)),[1 -.8],0,sfcast+4:sfcast+6);

cond.PIE_UTILS = obsdb.OBS_PIE_UTILS; %регулируемые тарифы (вернее, ЖКХ)
sspath = tseries(sfcast+25:efcast,m.PIE_UTILS_SS*ones(efcast-sfcast-25+1,1),get(cond.PIE_UTILS,'comment'));
obspath = resize(cond.PIE_UTILS,get(cond.PIE_UTILS,'first'):qq(2018,4));
cond.PIE_UTILS = smoothly_prolong(1000, 0.3, obspath, sspath);
p = exogenize(p, 'PIE_UTILS', sfcast:efcast);
p = endogenize(p, 'RES_PIE_UTILS', sfcast:efcast);

cond.OBS_PIE4 = tseries(sfcast, 100*log(1+0.01*(11.4)));% КСП по компонентам инфляции
p = exogenize(p, 'OBS_PIE4', sfcast);                   %выходим на 11.4 по 14 году
p = endogenize(p, 'RES_PIE_SERV_NO_UTILS', sfcast);     % за счет компонентыб которую досчитываем
% cond.OBS_LCPI_SERV = obsdb.OBS_LCPI_SERV; %услуги в целом (за счет нерег.)

cond.OBS_PIE = tseries(sfcast+1:sfcast+2,100*log(1+0.01*[28.0303 19.4697])); % 24.571 18.531   20.757 (14Q4) %%% 18.48 28.52 -> 19.88 28.83 ->23.74 53.14 39.07 -> 20.76 46.54 40.26
% p = exogenize(p, 'OBS_LCPI_SERV', sfcast:sfcast+1);
p = exogenize(p, 'OBS_PIE', sfcast+1:sfcast+2);
p = endogenize(p, 'RES_PIE_SERV_NO_UTILS', sfcast+1:sfcast+2);% КСП по компонентам инфляции

% cond.RES_PIE_SERV_NO_UTILS(sfcast+3:sfcast+3)=0.5*[-5.3774]; %ликвидация эффекта от курсовой динамики (накопившейся в 2014Q4-2015Q1cond.RES_PIE_SERV_NO_UTILS(sfcast+4:sfcast+6) = arf(tseries(sfcast+3,-3.4),[1 -.8],0,sfcast+4:sfcast+6);
% cond.RES_PIE_SERV_NO_UTILS(sfcast+4:sfcast+6) = arf(tseries(sfcast+3,cond.RES_PIE_SERV_NO_UTILS(sfcast+3)),[1 -.8],0,sfcast+4:sfcast+6);

%КСП по ВВП
cond.OBS_Y = obsdb.OBS_Y{sfcast}; 
cond.OBS_Y(sfcast:sfcast+2) = 100*log([1.793 1.732 1.726]);  %1.7926 1.7267 -> 1.7923 1.7078 -> 1.79 1.708 1.667->1.791 1.726 1.715
p = exogenize(p, 'OBS_Y', sfcast:sfcast+2);
p = endogenize(p, 'RES_Y_GAP', sfcast:sfcast+2);

% Внешний сектор
cond.OBS_PIEF = obsdb.OBS_PIEF{sfcast:end}; %Внешняя инфляция
p = exogenize(p, 'OBS_PIEF', sfcast:get(obsdb.OBS_PIEF,'end'));
p = endogenize(p, 'RES_PIEF', sfcast:get(obsdb.OBS_PIEF,'end'));

cond.OBS_Y_GAPF = obsdb.OBS_Y_GAPF{sfcast:end}; %Внешний "спрос"
p = exogenize(p, 'OBS_Y_GAPF', sfcast:get(obsdb.OBS_Y_GAPF,'end'));
p = endogenize(p, 'RES_Y_GAPF', sfcast:get(obsdb.OBS_Y_GAPF,'end'));

cond.OBS_RSF = obsdb.OBS_RSF; %Внешние ставки
sspath = tseries(sfcast+13:efcast,(m.RRF_SS+m.PIEF_SS)*ones(efcast-sfcast-13+1,1),get(cond.RSF,'comment'));
obspath = resize(cond.OBS_RSF,get(cond.OBS_RSF,'first'):get(obsdb.OBS_RSF,'end'));
cond.OBS_RSF = smoothly_prolong(1000, 0.3, obspath, sspath);
p = exogenize(p, 'OBS_RSF', sfcast:get(cond.OBS_RSF,'end'));
p = endogenize(p, 'RES_RRF_GAP', sfcast:get(cond.OBS_RSF,'end'));

% кросс-курс в перспективе - без изменений от последнего КСП
cond.USD_EUR = tseries(sfilt:sfcast+3,obsdb.USD_EUR(sfilt:sfcast+3));
cond.USD_EUR(sfcast+4:efcast) = cond.USD_EUR(sfcast+3)-.01;

%Цены на нефть
% cond.OBS_LPOIL = obsdb.OBS_LPOIL{sfcast:end};
cond.OBS_LPOIL(sfcast)=100*log((0.55 + 0.45./cond.USD_EUR(sfcast))*75.9);%78
% cond.OBS_LPOIL(sfcast+1:sfcast+4)=100*log((0.55 + 0.45./cond.USD_EUR(sfcast+1:sfcast+4))*40);
% cond.OBS_LPOIL(sfcast+5)=100*log((0.55 + 0.45./cond.USD_EUR(sfcast+5))*50);
% cond.OBS_LPOIL(sfcast+6:sfcast+12)=100*log((0.55 + 0.45./cond.USD_EUR(sfcast+6:sfcast+12))*60);
cond.OBS_LPOIL(sfcast+1:sfcast+4)=100*log((0.55 + 0.45./cond.USD_EUR(sfcast+1:sfcast+4)).*[50 54 60 64]');
cond.OBS_LPOIL(sfcast+5:sfcast+8)=100*log((0.55 + 0.45./cond.USD_EUR(sfcast+5:sfcast+8)).*(linspace(64,69,4)'));
% cond.OBS_LPOIL(sfcast+13:sfcast+14)= arf(cond.OBS_LPOIL{sfcast+12},[1 -1.025],0,sfcast+13:sfcast+14);
p = exogenize(p, 'OBS_LPOIL', sfcast:get(cond.OBS_LPOIL,'end'));
p = endogenize(p, 'RES_LTOT_GAP', sfcast:get(cond.OBS_LPOIL,'end'));

cond.TUNE_LTOT_GAP(sfcast+13:sfcast+20) = 0;
% cond.TUNE_LTOT_GAP(qq(2019,1):sfcast+20) = 0;
p = exogenize(p, 'TUNE_LTOT_GAP', qq(2018,1):sfcast+20);
p = endogenize(p, 'RES_LTOT_GAP', qq(2018,1):sfcast+20);


%% сценарные предпосылки прогноза - суждения

% цель по инфляции


% cond.TUNE_PIE_TAR(sfcast:sfcast+8) =  100*log([5.0...
%                                            0+[4.875 4.75 4.625 4.5]...
%                                            0+[4.375 4.25 4.125 4.0]]/100+1);
%       
cond.TUNE_PIE_TAR(sfcast:sfcast+12) =  100*log(linspace(round(5.0+(8.0661+10.7547+5.3774)/4),4,13)/100+1);                                     
p = exogenize(p, 'TUNE_PIE_TAR', sfcast:get(cond.TUNE_PIE_TAR,'end'));
p = endogenize(p, 'RES_PIE_TAR', sfcast:get(cond.TUNE_PIE_TAR,'end'));

% околонулевой разрыв по курсу                                     
% cond.RES_DOT_LZ_EQ(sfcast)=2;

% cond.TUNE_LZ_GAP(sfcast) = 0;
% p = exogenize(p, 'TUNE_LZ_GAP', sfcast);
% p = endogenize(p, 'RES_DOT_LZ_EQ', sfcast:sfcast+2);

% cond.TUNE_LZ_GAP(sfcast+2) = 24.69/2;
% cond.TUNE_LZ_GAP(sfcast+3) = 24.69/4;
% cond.TUNE_LZ_GAP(sfcast+4) = 0;
% p = exogenize(p, 'TUNE_LZ_GAP', sfcast+2:sfcast+4);
% p = endogenize(p, 'RES_LS', sfcast+2:sfcast+4);
% cond.RES_DOT_LZ_EQ(sfcast+5:sfcast+9) = arf(tseries(sfcast+4, -5.15),[1 -.4],0,sfcast+5:sfcast+9);

% прекращение улучшения условий торговли
cond.DOT_LTOT_EQ(sfcast:sfcast+17) = fdb2.mean.DOT_LTOT_EQ(sfcast:sfcast+17);
% cond.DOT_LTOT_EQ(sfcast:sfcast+17) = filtdb.mean.DOT_LTOT_EQ(sfcast:sfcast+17);
p = exogenize(p, 'DOT_LTOT_EQ', sfcast:get(cond.DOT_LTOT_EQ,'end'));
p = endogenize(p, 'RES_DOT_LTOT_EQ', sfcast:get(cond.DOT_LTOT_EQ,'end'));

% страновая премия
% obsdb.OBS_DOT_PREM = 100*log(1+.01*tmpdb.CDS)-100*log(1+.01*tmpdb.CDS{-1});
prem=[2.266 3.267 6.200 7.000]; %
prem=100*log(1+.01*prem(2:end))-100*log(1+.01*prem(1:end-1));
cond.OBS_DOT_PREM(sfcast:sfcast+2) = obsdb.OBS_DOT_PREM(sfcast:sfcast+2); % ниже уровень потом
cond.OBS_DOT_PREM(sfcast:sfcast+2)=prem;
% cond.OBS_DOT_PREM(sfcast:sfcast+12)=prem;
p = exogenize(p, 'OBS_DOT_PREM', sfcast:get(cond.OBS_DOT_PREM,'end'));
p = endogenize(p, 'RES_PREM', sfcast:get(cond.OBS_DOT_PREM,'end'));

% шок курса "вдобавок" к страновой премии (для объяснения динамики курса)
cond.OBS_LS(sfcast) = obsdb.OBS_LS;
cond.OBS_LS(sfcast:sfcast+1) = 100*log([52.87 68.5]); %52.93 80
p = exogenize(p, 'OBS_LS', sfcast:get(cond.OBS_LS,'end'));
p = endogenize(p, 'RES_LS', sfcast:get(cond.OBS_LS,'end'));

% % % % % % cond.RES_LS(sfcast) = tseries(sfcast+2:sfcast+4,0);  %интервенции
% % % % % % cond.RES_LS(sfcast+2) = -3.360124; %
% % % % % % cond.RES_LS(sfcast+3) = 4.815792; %
% % % % % % cond.RES_LS(sfcast+4) = -18.027386; %
% % % % % % p = exogenize(p, 'RES_LS', sfcast+2:sfcast+4);
% % % % % % p = endogenize(p, 'OBS_LS', sfcast+2:sfcast+4);

% cond.RES_LS(sfcast:sfcast+12)=[2.6 2.2 1.9 1.6 1.4 1.2 1 .9 .7 .6 .5 .4 .4];% как в сентябре
% cond.std_RES_LS(sfcast+2:sfcast+5)=0;

% % % % p = endogenize(p, 'RES_DOT_LZ_EQ', sfcast); %шамаааанство)))  :sfcast+1
% % % % p = endogenize(p, 'RES_LS', sfcast+1);


% cond.RES_PREM = tseries(sfcast:sfcast+2,1.5*[2 1.3 1]+1i*[0 -1 -0.5]); %[5i 1.8i 1.4i]
% cond.RES_PREM(get(cond.RES_PREM,'end')+1:sfcast+11)=arf(tseries(get(cond.RES_PREM,'end'),cond.RES_PREM(get(cond.RES_PREM,'end'))),[1 -.5],0,get(cond.RES_PREM,'end')+1:sfcast+11);
% cond.std_RES_PREM = tseries(sfcast:sfcast+11,0);

% снижение потенциального ВВП - бодрее, чем в феврале
cond.RES_DOT_Y_EQ(sfcast:sfcast+4)=.8*[-.2 -.5 -.5 -.3 -.2];
cond.std_RES_DOT_Y_EQ = tseries(sfcast:sfcast+4,0); % бодрее спад вначале

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % cond.TUNE_DOT_Y_EQ = tseries(sfcast+6:sfcast+18,1.25); % ниже уровень потом
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % p = exogenize(p, 'TUNE_DOT_Y_EQ', sfcast+6:sfcast+18);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % p = endogenize(p, 'RES_DOT_Y_EQ', sfcast+6:sfcast+18);

cond.TUNE_DOT_Y_EQ = forecastBASELINE.mean.DOT_Y_EQ; % ниже уровень потом
p = exogenize(p, 'TUNE_DOT_Y_EQ', sfcast:get(cond.TUNE_DOT_Y_EQ,'end'));
p = endogenize(p, 'RES_DOT_Y_EQ', sfcast:get(cond.TUNE_DOT_Y_EQ,'end'));

cond.TUNE_DOT_LZ_EQ(sfcast:sfcast+4) = [[50 30 15]*0.73 3 2]; %[50 30 15]*.8 3 2
p = exogenize(p, 'TUNE_DOT_LZ_EQ', sfcast:get(cond.TUNE_DOT_LZ_EQ,'end'));
p = endogenize(p, 'RES_DOT_LZ_EQ', sfcast:get(cond.TUNE_DOT_LZ_EQ,'end'));
% % % cond.TUNE_LZ_GAP(sfcast+2) = 17.35/2;
% % % cond.TUNE_LZ_GAP(sfcast+3) = 17.35/4;
% % % cond.TUNE_LZ_GAP(sfcast+4) = 0;
% % % p = exogenize(p, 'TUNE_LZ_GAP', sfcast+2:sfcast+4);
% % % p = endogenize(p, 'RES_LS', sfcast+2:sfcast+4);

% экзогенная ставка до 2014Q4
cond.OBS_RS(sfcast) = 100*log(1+.01*[11.392]); %11.959
% cond.OBS_RS(sfcast:sfcast+0) = obsdb.OBS_RS;
% cond.OBS_RS(sfcast+1:sfcast+4) = 100*log(1+.01*12);
p = exogenize(p, 'OBS_RS', sfcast:sfcast+0);
p = endogenize(p, 'RES_RS', sfcast:sfcast+0);

cond.OBS_RS_MARKET(sfcast:sfcast+0) = 100*log(1+.01*[12.6]); %13.17
p = exogenize(p, 'OBS_RS_MARKET', sfcast:sfcast+0);
p = endogenize(p, 'RES_RS_MARKET', sfcast:sfcast+0);

%% сценарные предпосылки прогноза - техническое и разглаживание

% зануление шоков в уравнениях сборки
cond.std_RES_LCPI= tseries(sfcast:sfcast+3,0);
cond.std_RES_LCPI_SERV= tseries(sfcast:sfcast+11,0);
    
% сглаживание шока ВВП - ПОДОБРАТЬ КОНСТАНТУ!!!
% Это наши предположения об эффектах от расходования фондов:
% cond.RES_Y_GAP(sfcast+2:sfcast+12) = arf(tseries(sfcast+1,-1.44),[1 -.8],0,sfcast+2:sfcast+12);
cond.RES_Y_GAP(sfcast+3:sfcast+12) = arf(tseries(sfcast+2,-0.44*0.8),[1 -.8],0,sfcast+3:sfcast+12)... эта строчка - про эффект санкций
                                    +arf(tseries(sfcast+2,1.17),[1 -.8],0,sfcast+3:sfcast+12);        %эта строчка - про нелинейность зависимости от нефти
% cond.RES_Y_GAP(sfcast+1) = -1;
cond.RES_Y_GAP_CRISIS(sfcast+1:sfcast+12) = 0.5*arf(tseries(sfcast,1.6*1i),[1 -.83],0,sfcast+1:sfcast+12);
% cond.RES_Y_GAP(qq(2015,2):qq(2015,3)) = cond.RES_Y_GAP(qq(2015,2):qq(2015,3));

% сглаживание шока курса (кагбе-премии)
% cond.RES_LS(sfcast+1:sfcast+12) = arf(tseries(sfcast,5.5),[1 -.55],0,sfcast+1:sfcast+12);
% сглаживание шока рыночной ставки - ПОДОБРАТЬ КОНСТАНТУ!!!
% cond.RES_RS_MARKET(sfcast:sfcast+11)=arf(tseries(sfcast-1,.12),[1 -.75],0,sfcast:sfcast+11);
% cond.std_RES_RS_MARKET = tseries(sfcast:sfcast+11,0);

%% базовый прогноз
    % сохранение сценарных предпосылок и КСП
    save([out_dir 'bl_cond' ntfsfx '.mat'],'cond');
    
    % прогноз
     simdb = simulate(m, cond, sfcast:efcast, 'plan', p, 'anticipate', true);
     fdb.mean = dbextend(initdb.mean,simdb); 
     fdb.mean = reporting(m, fdb.mean, sfilt:efcast);
   
% дополнительные переменные
     fdb.mean.RPT_LPOIL_USD_F = comment(exp(0.01*fdb.mean.LPOIL)/(.55+.45/cond.USD_EUR),...
        'Цена на нефть,\it~долл~за~брл');
     fdb.mean.RPT_LTOT = comment(dbeval(fdb.mean,...
        'exp(0.01*(LTOT+LCPIF(qq(2005,1))-LPOIL(qq(2005,1))))'),...
        'Реальная цена на нефть,\it~2005q1=1');
     fdb.mean.RPT_LTOT_EQ = comment(dbeval(fdb.mean,...
        'exp(0.01*(LTOT_EQ+LCPIF(qq(2005,1))-LPOIL(qq(2005,1))))'),...
        'Реал. равновес. цена на нефть,\it~2005q1=1');

 % сохранение результатов
    save([out_dir 'bl_fcast' ntfsfx '.mat'],'fdb');
    dbsave([out_dir 'bl_fcast' ntfsfx '.csv'],fdb.mean,get(fdb.mean.PIE,'range'),...
        'format','%f','class',false);

    
    % сохранение модели, используемой для прогноза
    save([out_dir 'model_bl_fcast.mat'],'m');     
     return;

%% напоминалка и вывод отдельных результатов в командном окне 

% создание базы данных для годного к отображению периода
my = dbclip(fdb.mean,sfcast-1:sfcast+11);

% выбор переменных
% [my.RPT_PREM my.RPT_LS my.RPT_RS my.RPT_PIE4 my.RPT_DOT_Y my.RPT_DOT_Y4 my.RPT_Y_GAP]

% напоминалко
fprintf('\nВсе результаты лежат в папке <a href="%s">%s</a>\n',...
    [cd '\' out_dir],out_dir);
