%% Start IRIS if it is not already running:
a1_init;

%% Work folders:
res_dir = ['results' filesep 'filter_data'];

%% Load the model:
m = model('reforma.mod','linear',true);

%% Changes to the model parameters can be put here

%% Solve the model
m = solve(m);
m = sstate(m);
% return

%% Load observed data:
load(fullfile(res_dir,'processed_data.mat'),'db');

% obsdb_hard = struct(); % obsdb without tunes, hard data only
% obsdb      = struct(); % obsdb with all the tunes

% add empty tseries for missing tunes (to avoid warnings)
yList2add = get(m,'yList');
for flx = yList2add
  obsdb.(flx{1}) = tseries();
end



obsdb.obs_y = db.gdp;
obsdb.obs_pi = 4*diff(db.cpi_mom_sa);
obsdb.obs_rs =db.rs;
obsdb.obs_s = db.s;
obsdb.obs_prem = db.prem;
obsdb.obs_oil = db.oil;
obsdb.obs_y_gap_f= hpf2(db.gdp_f);
 obsdb.obs_rs_f = db.rs_f;
 obsdb.obs_pi_f  = 4*diff(db.cpi_f_mom_sa);
 
 %% Assumptions
% Inflation target:
obsdb.tune_pi_tar(qq(2015,3):qq(2018,3))= [100*log(5.5*0.01+1)... 
                                                  100*log(5.375*0.01+1) 100*log(5.25*0.01+1) 100*log(5.125*0.01+1) 100*log(5.0*0.01+1)...
                                                  100*log(4.875*0.01+1) 100*log(4.75*0.01+1) 100*log(4.625*0.01+1) 100*log(4.5*0.01+1)...
                                                  100*log(4.375*0.01+1) 100*log(4.25*0.01+1) 100*log(4.125*0.01+1) 100*log(4.0*0.01+1)];
%                                               
% obsdb.tune_tot_gap(qq(2013,2))= 0;
% obsdb.tune_tot_gap(qq(2013,3))= 0;
% obsdb.tune_tot_gap(qq(2015,2))= 0;
% Output gap:
obsdb.tune_y_gap(qq(2013,3):qq(2013,4)) = 0;
obsdb.tune_y_gap(qq(2015,3):qq(2015,4)) = -1.5;     
                                              
save(fullfile(res_dir,'obs_data.mat'),'obsdb');

%% Quasi-parametrisation, std multipliers:
stds     = get(m,'std');
list_std = fieldnames(stds);
mult     = dbbatch(stds,'$0','tseries(sfilt:efcast,1)','namelist',list_std); 

for ii = 1:numel(list_std)
   if ~isempty(strfind(list_std{ii},'_obs_'))
      obs_var = strrep(list_std{ii},'std_eps_obs_','');
      mult.(list_std{ii})(sfilt:efcast) = 0;
   end
end

% Higher vol of target before crisis
mult.std_eps_pi_tar(sfilt:qq(2009,2)) = 2;

std_tune = dbbatch(mult,'$0','stds.$0*mult.$0','namelist',list_std);  

%% Kalman filter
if length(m)==1
   P = get(m,'param');
   m = assign(m,P);
   [mfilt,filtdb,se2,delta,pe] = filter(m,obsdb,sfilt:efcast,'std',std_tune);
%    report_likelihood_contrib(filtdb.mean,std_tune,sfilt:efcast,list_std,se2);
else
   [mfilt,filtdb,se2,delta,pe] = filter(m,obsdb,sfilt:efcast);
end

save(fullfile(res_dir,'model_filtration.mat'), 'mfilt');
% save(fullfile(res_dir,'std_tune.mat'), 'std_tune');

%% Save results

save(fullfile(res_dir,'filtered_dbase.mat'), 'filtdb')
my = dbclip(filtdb.mean,sfcast-8:sfcast+7);

%% 
fprintf('\nAll results are <a href="%s">%s</a>\n', [cd '\' res_dir],res_dir);
return