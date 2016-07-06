%%%%%%%%%%%%%% reforma %%%%%%%%%%%%%%

%%  ver. 0.7 as of 2016-01-14
% Add into other key equations these blocks:
% * ToT block        V
% * FX block         

%%  ver. 0.6 as of 2016-01-14
% Added terms of trade block. 

%%  ver. 0.5 as of 2016-01-13
% Deleted multiple Phillips curves from the model, now left just one for
% simplicity.

%% ver. 0.4 as of 2016-01-13
% iav: 
% * ABorodin tweaked the UIP so that it should now be working. 

%% ver. 0.3 as of 2016-01-08
% iav: 
% * Terms of Trade             
% * FX                                  ?

%% ver. 0.2 as of 2016-01-08
% iav: 
% * IS curve                        V
% * Phillips curve               V
% * Interest rates               V
% * FX                                  ?
% * Foreign economy        V

%% ver. 0.1 as of 2016-01-07
% iav: I've typed in all the equations but the model is not functional. The problem is that the solution is singular.
%		Should check the equations one by one..
%% TODO:
% * GDP decomposition
% * labour market 
% * regulators and general public expectations of inflation
% * relative prices gaps
% * wedge between market money market rate and policy key rate

%% ver. 0.0 as of 2016-01-06
% This is a verbatum copy of the model in 
% Dizioli, Schmittmann (2015) A Macro-Model Approach to Monetary Policy Analysis and Forecasting for Vietnam, IMF, WP/15/273
% The difference is
% * inflation is devided into 4 components: food, non-food, non-regulated services, regulated services
% * the foreign block is not forecasted but is fully exogenous

%%%%%%%% Notation %%%%%%%%
%% _qoq - quarterly change
%% _yoy - yearly change 
%% e_ - expectations
%% obs_ - observed
%% _f - foreign variable


!transition_variables
%% Output growth
	'GDP, const prices' y 
	'GDP trend' y_eq
	'Output gap' y_gap
	'GDP, % QoQ' y_qoq
	'GDP, % YoY' y_yoy
	'GDP trend, % QoQ' y_eq_qoq
% 	'GDP, rolling 4Q YoY, %' y_yoy_ann
	'Expected 1Q ahead output gap' e_y_gap
% Inflation
	'Level of CPI' cpi
	'CPI inflation' pi
	'CPI inflation, YoY' pi_yoy
	'CBRs inflation target' pi_tar
	'' e_pi
	'' e_pi_tar
    e_pi_dev
    
% Interest rates and monetary policy	
	'Key rate' rs
	'Real key rate' rr
	'Real key rate, equillibrium' rr_eq
	'Real rate gap' rr_gap
	'Country risk premium' prem

% Exchange rates
	'Nominal exchange rate' s
	'Nominal exchange rate, % QoQ' s_qoq
	'Real exchange rate' z
	'RER, equillibrium' z_eq
    z_gap
	'RER, % QoQ'  z_qoq
	'RER, equillibrium' z_eq_qoq
	'Expectations of RER change' e_z
	e_z_eq_qoq
    'Exchange rate expectations' e_s
    'Exchange rate expectations' ew_s

% Terms of trade block
	'Urals price, USD' oil
	'Terms of trade' tot
	'Terms of trade, trade' tot_eq
	'Terms of trade, gap' tot_gap
	'Terms of trade, trend, qoq' tot_eq_qoq
    
% Foreign block
	'Foreign real interest rate' rr_f
	'Foreign real interest rate, equillibrium' rr_f_eq
	'Foreign real interest rate gap' rr_f_gap
	'Foreign nominal key interest rate' rs_f
	'Foreign output gap' y_gap_f
	'Foreign inflation' pi_f
    cpi_f
	e_pi_f
	
%% GDP. Output gap and equillibrium output growth:
!transition_equations
	y_gap = a_y_gap_1*e_y_gap + a_y_gap_2*y_gap{-1}-a_y_gap_3*rr_gap{-1}+a_y_gap_4*z_gap{-1}+a_y_gap_5*y_gap_f{-1}+a_y_gap_6*tot_gap+eps_y_gap;
    y_eq_qoq = a_y_eq_qoq1*y_eq_qoq{-1} + (1-a_y_eq_qoq1)*(-a_y_eq_qoq2*z_eq_qoq+y_qoq_ss) + eps_y_eq_qoq;
	
	e_y_gap=y_gap{+1};
	y=y_eq+y_gap;
	y_qoq=4*(y-y{-1});
	y_yoy=(y-y{-4});
    y_eq = y_eq{-1} + 0.25*y_eq_qoq;
	
	!transition_shocks
	'IS shock' eps_y_gap
	'Equillibrium growth rate shock' eps_y_eq_qoq
	!parameters
	a_y_gap_1=0.7
 	a_y_gap_2=0.2
	a_y_gap_3=0.3
	a_y_gap_4=0.1
	a_y_gap_5=0.6
    a_y_gap_6=0.1
	a_y_eq_qoq1=0.2
    a_y_eq_qoq2=0.3
	y_qoq_ss=100*log(1.5*0.01+1)

% %% Phillips curves:
	!transition_equations
%   We delete variable related to the foreign block, ie exchange rate:
    pi=a_pi_1*e_pi+(1-a_pi_1)*pi{-1}+a_pi_3*y_gap{-1}+a_pi_4*(z_gap-z_gap{-1})+eps_pi;
%  Expectations (HERE I SHOULD PUT s_qoq)
%     e_pi=pi_yoy{+4};
% 	e_pi=a_e_pi_1*pi_yoy{+1}+a_e_pi_2*pi{-1}+(1-a_e_pi_1-a_e_pi_2)*s_qoq{-1};
	e_pi=a_e_pi_1*pi_yoy{+1}+(1-a_e_pi_1)*pi{-1};
 	%  Time-varying target:
	pi_tar=a_pi_tar*pi_tar{-1}+(1-a_pi_tar)*pi_tar_ss+eps_pi_tar;
	
    cpi=cpi{-1}+0.25*pi;
	pi_yoy=(cpi - cpi{-4});
	!transition_shocks
	'Inflation shock' eps_pi
    'Inflation target shock' eps_pi_tar
	!parameters
	pi_tar_ss=100*log(6.0*0.01+1)
	pi_rs_ss=100*log(0.7*pi_tar_ss*0.01+1)
    a_pi_tar=0.6
	% Non-food parameters
	a_pi_1=0.6
	a_pi_3=0.3
    a_pi_4=0.3
%  Expectations   
    a_e_pi_1=0.6;
	a_e_pi_2=0.3;

% %% Monetary policy
!transition_equations
% %  Monetary policy rule (this from the paper, not from the CBR's MPR)
% 	rs=a_rs_1*rs{-1}+(1-a_rs_1)*(rr_eq+pi_yoy+a_rs_2*(pi_yoy{+4}-pi_tar{+4})+a_rs_3*y_gap)+a_rs_4*(z-z{-1})+eps_rs;
%  Without the real rate gap:
    rs=a_rs_1*rs{-1}+(1-a_rs_1)*(rr_eq+e_pi+a_rs_2*e_pi_dev+a_rs_3*y_gap)+eps_rs;
% %  Neutral real interest rate 2.5 should actually be prem, and z_qoq_eq which should be dependant on the terms of trade shifts should also be here:
    rr_eq=a_rr_eq*rr_eq{-1}+(1-a_rr_eq)*(prem+z_eq_qoq)+eps_rr_eq;

    rr=rs-e_pi;
	rr_gap=rr-rr_eq;
    e_pi_tar=pi_tar{+4};
    e_pi_dev=(e_pi-e_pi_tar);
    e_z_eq_qoq=z_eq_qoq{+1};
	!transition_shocks
	'Monetary policy shock' eps_rs
	'Neutral interest rate shock' eps_rr_eq
% 	
	!parameters
	a_rs_1=0.4
	a_rs_2= 2.5
	a_rs_3=2.5
    % a_rs_4=0.2 % real exchange rate coefficient
	a_rr_eq=0.5
 	
    % % % % % % % % % % % % % % % % % % % % % % % %% UIP
!transition_equations
% This is the QPM version of the UIP:
    s = ew_s - (rs-rs_f-prem)/4+eps_s/4;
    ew_s=a_e_s*e_s+(1-a_e_s)*(s{-1}-(pi-pi_f+e_pi-e_pi_f-2*z_eq_qoq)/4);  
  
    e_s=s{+1};
    e_z=z{+1}; % 
    % % % % % % % % % % % % % % % % % % % % % %     % Change in the equillibrium z, should add here a section on ToT:
   z_eq_qoq    = a_z_eq_qoq1*z_eq_qoq{-1} + (1-a_z_eq_qoq1)*(-a_z_eq_qoq2*tot_eq_qoq) + eps_z_eq_qoq;    
    % % % % % % % % % % % % % % % % % % % % % %      % Country risk premium 
   prem=a_prem*prem{-1} + (1-a_prem)*prem_ss + eps_prem;
     % % % % % % % % % % % % % % % % % % % % % %     e_s=s{+1};
    s=s{-1}+0.25*s_qoq;
    s_qoq= z_qoq +pi-pi_f;
    z=z{-1}+0.25*z_qoq;
	z_eq=z_eq{-1} + 0.25*z_eq_qoq;
    z_gap = z-z_eq;
    % % % % % % % % % % % % % % % % % % % % % % 
	!transition_shocks
	eps_s
	eps_prem
	eps_z_eq_qoq
	!parameters
	a_e_s=0.1
	prem_ss=2.0
	a_prem=0.2
	a_z_eq_qoq1=0.3
	a_z_eq_qoq2=0.7
	tot_qoq_ss=6.0
    
    %% Terms of trade
    !transition_equations
    tot_eq_qoq = a_tot_eq_qoq*tot_eq_qoq{-1}+(1-a_tot_eq_qoq)*tot_qoq_ss+eps_tot_eq_qoq;
    tot_gap     = a_tot_gap*tot_gap{-1} + eps_tot_gap;
    
    oil=tot+cpi_f;
    tot=tot_eq+tot_gap;
    tot_eq=tot_eq{-1}+0.25*tot_eq_qoq;
    
    !transition_shocks
    eps_tot_eq_qoq
    eps_tot_gap
    !parameters
    a_tot_eq_qoq=0.7;
    a_tot_gap=0.6
    
%% Foreign economy
!transition_equations
	y_gap_f=a_y_gap_f*y_gap_f{-1}+0.0+eps_y_gap_f;
	pi_f=a_pi_f*pi_f{-1}+(1-a_pi_f)*pi_f_ss+eps_pi_f;
	rr_f_eq=a_rr_f_eq*rr_f_eq{-1}+(1-a_rr_f_eq)*rr_f_ss+eps_rr_f_eq;
	rr_f_gap=a_rr_f_gap*rr_f_gap{-1}+eps_rr_f_gap;
    
    cpi_f=cpi_f{-1}+0.25*pi_f;
% 	
    rs_f=rr_f+e_pi_f;
	e_pi_f=pi_f{+1};
	rr_f=rr_f_eq+rr_f_gap;
	!transition_shocks
	eps_y_gap_f
	eps_pi_f
	eps_rr_f_eq
	eps_rr_f_gap
	!parameters
	a_y_gap_f=0.8
	a_pi_f=0.8
	pi_f_ss=100*log(2.0*0.01+1)
	a_rr_f_eq=0.4
	rr_f_ss=0.0
	a_rr_f_gap=0.1


%% Measurement part
!measurement_variables
	'Observed GDP' obs_y
%    
	'Observed inflation' obs_pi
% 	
	'Observed key rate' obs_rs
% 	
	'Observed nominal exchange rate' obs_s
	'Country risk premium' obs_prem
% 
    'Observed Urals price' obs_oil
% 	
	'Observed foreign GDP gap' obs_y_gap_f
	'Observed foreign nominal key rate' obs_rs_f
	'Observed foreign inflation rate' obs_pi_f
	
!measurement_equations
	obs_y=y;
% 
    obs_pi=pi;
% 
	obs_rs=rs;
% 
	obs_s=s;
	obs_prem=prem;
%
    obs_oil = oil;
% 
	obs_y_gap_f=y_gap_f;
	obs_rs_f=rs_f;
	obs_pi_f=pi_f;

% Reporting equations
% !reporting_equations
% !import(./tunes.mod)
% !import(./reportingvars.mod)