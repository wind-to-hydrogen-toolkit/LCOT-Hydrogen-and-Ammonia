clear all
close all
clc
clf

fig1 = figure('numbertitle','off',...
    'units','normalized',...
    'menubar','none',...
    'position',[0.2 0.1 0.6 0.85],...
    'color',[1 1 1],...
    'visible','off');

%Inputs
Power = [200:100:3000];                      %Electrolyzer plant size [MW] 1443
Dist = 100:50:500;                      %Distance [km]
L = 0:1:30;                               %Overall plant lifetime [years]
dr = 0.08;                                %Discount rate = 8%
Elec_prod = 0.0055;                       %Electrolyzer production rate [kg/s/MW] Siemens(Silyzer300)
CF = 1;                                   %Electrolyzer capacity factor = 100%
Urate_p = 0.75;                           %Utilization rate of the pipeline = 75% (IEA report)
OPEX_rate_p = 0.02;                       %OPEX = 2% of CAPEX
OPEX_rate_liq = 0.04;                     %OPEX = 4% of CAPEX
OPEX_rate_stor = 0.04;                    %OPEX = 4% of CAPEX
OPEX_rate_ship = 0.04;                    %OPEX = 4% of CAPEX
OPEX_rate_HVDC = 0.04;                    %OPEX = 4% of CAPEX
V_ship = 30;                              %Velocity of ship (km/h)
Storage_cost = 85;                        %Storage cost (EUR/kgH2)

%---------------------------------------------------------------------------------------------------------------------
figure(1)
break_dist = [];
for k = Power
    LCOH2_p_m = [];                       %Creating matrix

    LCOH2_sto_m = [];                     %Creating matrix
    LCOH2_ship_m = [];                    %Creating matrix

    for i = Dist
        ratio_Ammonia = 17/3;
        %Control variables
        Prod_h = k*Elec_prod*CF*3600;      %Amount of hydrogen produced per hour [kg/h]
        Prod_d = Prod_h*24;                %Amount of hydrogen produced per day [kg/d]
        Prod_y = Prod_d*350;               %Amount of hydrogen produced per year [kg/y]
        Prod_l = Prod_y*max(L);            %Amount of hydrogen produced over the plant lifetime [kg]

        Prod_A_h = Prod_h*ratio_Ammonia;
        Prod_A_d = Prod_d*ratio_Ammonia;
        Prod_A_y = Prod_y*ratio_Ammonia;
        Levelized_prod = 0;                %Initiating variable
        Levelized_prod_pipe = 0;

        for j = L
            Levelized_prod = Levelized_prod + (Prod_A_y/((1+dr)^j));  %NPV of production
        end
        for j = 1:1:40
            Levelized_prod_pipe = Levelized_prod_pipe + (Prod_A_y/((1+dr)^j));  %NPV of production
        end
%======================================= PIPELINES =======================================%
        D = power(Prod_A_h/1000/15.245,1/1.6468)*5 + 5;
        n_pump = floor(i/128);
        CAPEX_p_km = D*857000*0.9/25; %Capital cost per km (EUR/km]
        CAPEX_com = (n_pump*2000000)*0.9;
        CAPEX_p = 2*((i*CAPEX_p_km) + CAPEX_com); % %CAPEX multiplied by distance brought to present value 2019 (source 2010)

%======================================== SHIP ========================================%

        %-------------------- CAPEX Storage  --------------------%
        num_ship = 1;
        Cap_ship = 53000000;
        Capex_tank = 68*10^6*0.9;
        Cap_tank = 34100000; % kg NH3
        timetravel_pership = 2*i/V_ship;
        check = timetravel_pership*Prod_A_h;
        if check > Cap_ship
          num_ship = ceil(check/Cap_ship)
        endif

        Num_tank = ceil(Cap_ship*num_ship/Cap_tank);
        CAPEX_stor = Num_tank*Capex_tank;

        %--------------------  Charter rates  --------------------%

        CAPEX_ship = 85*10^6*num_ship*0.9;

        Price_HFO = 571.39*0.9;  % Heavy fuel oil 571.39 usd/ton  https://www.iea.org/reports/key-world-energy-statistics-2020/prices
        E_density_HFO = 39;      % MJ/kg https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html
        E_cost_unit = Price_HFO/1000/39; % Euro/MJ
        Fuel_ship_km = E_cost_unit*2500; % Fuel use 2500 MJ/km
        Fuel_cost = 2*i*Fuel_ship_km*ceil(Prod_A_y/Cap_ship); % 30euro/km
        %CAPEX_ship = (Price_Bench*(Ship_cap/(Cap_Bench*1000))^0.623)*1000000;
        %Charter_price_LH2 = CAPEX_ship*Ship_cap_freight_ratio;


        %---Calculating the NPV of OPEX for pipeline and ship costs -------------------------------------%

        OPEX_p = 0;            %Initiating and restarting value for calculating NPV for next distance
        OPEX_ship = 0;         %Initiating and restarting value for calculating NPV for next distance


        for j = 1:1:40
          OPEX_p = OPEX_p + (CAPEX_p*OPEX_rate_p + (n_pump*408000*0.9))/((1+dr)^j);               %NPV of OPEX for pipeline
        endfor

        for j = L

            OPEX_ship = OPEX_ship + (CAPEX_stor*OPEX_rate_stor+...            %OPEX storage
                                     CAPEX_ship*OPEX_rate_ship+...            %OPEX tanker
                                     Fuel_cost)/((1+dr)^j);           %Charter rates

        end

        LCOH2_p = (CAPEX_p+OPEX_p)/Levelized_prod_pipe;           %Levelized cost of hydrogen for pipeline transport[EUR/kgH2]
        LCOH2_p_m = [LCOH2_p_m,LCOH2_p];

        LCOH2_ship = (CAPEX_stor+OPEX_ship+CAPEX_ship)/Levelized_prod;  %Levelized cost of hydrogen for ship transport (storage + ship + freight) [EUR/kgH2]
        LCOH2_ship_m = [LCOH2_ship_m,LCOH2_ship];  %Levelized cost of hydrogen for ship transport (Liquefaction + storage + ship + freight) [EUR/kgH2]

    end
    tmp = abs(LCOH2_p_m-LCOH2_ship_m);
    idx = find(tmp==min(tmp));
    break_dist = [break_dist, Dist(idx)];

    hold on
##    plot (Dist,LCOH2_p_m,'Color','k','linewidth',1.5)
##    plot (Dist,LCOH2_ship_m,'--','Color','b','linewidth',1.5)
##    %plot (Dist_HVDC,LCOH2_HVDC_m,'Color','r')
##    %plot (Dist,LCOH2_liq_m,'Color',[0, 0.5, 0],'linewidth',1.5)
##    text(max(Dist)-25,max(LCOH2_p_m)+0.001,[num2str(k),' MW'],'fontsize',[18],...
##    'fontname','times','color',[0 0 0],'HorizontalAlignment','right')
##    text(max(Dist)-25,max(LCOH2_ship_m)+0.005,[num2str(k),' MW'],'fontsize',[18],...
##    'fontname','times','color','b','HorizontalAlignment','right')
end
##LCOT_p = LCOH2_p_m(41)
##LCOT_ship = LCOH2_ship_m(41)
##set (gca,'fontsize',20,'fontname','times','box','on','xlim',[100 max(Dist)])
##xlabel('Transmission distance (km)','fontname','times')
##ylabel('LCOT (EUR/kgNH_{3})','fontname','times')
##legend ('Pipeline','Ammonia tankers'...
##        ,'fontsize',[20],'Location','Northwest')

##figure(2)
plot (break_dist,Power,'r','linewidth',1.5)
##
xlabel('Break-Even distance (km)','fontname','times','fontsize',20)
ylabel('Electrolyzer capacity (MW)','fontname','times','fontsize',20)
set (gca,'fontsize',20,'fontname','times','box','on','xlim',[100 max(break_dist)])
%set (gca,'fontsize',20,'fontname','times','box','on','xlim',[200 max(Power)])
t = cputime

