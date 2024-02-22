clf
%Inputs
Power = 500:100:3500;                     %Electrolyzer plant size [MW]
Dist = 100:50:1000;              %Distance [km]
ZZ = zeros (numel (Power), numel (Dist));
L = 0:1:30;                               %Overall plant lifetime [years]
dr = 0.08;                                %Discount rate = 8%
Elec_prod = 0.0055;                       %Electrolyzer production rate [kg/s/MW] Siemens(Silyzer300)
CF = 1;                                   %Electrolyzer capacity factor = 100%
Urate_p = 0.75;                           %Utilization rate of the pipeline = 75% (IEA report)
OPEX_rate_p = 0.02;                       %OPEX = 2% of CAPEX
OPEX_rate_liq = 0.04;                     %OPEX = 4% of CAPEX
OPEX_rate_stor = 0.04;                    %OPEX = 4% of CAPEX
OPEX_rate_ship = 0.04;                    %OPEX = 4% of CAPEX

V_ship = 30;                              %Velocity of ship (km/h)
Storage_cost = 85;                        %Storage cost (EUR/kgH2)

%Relation between charter rates and LH2 vessel capacity according to benchmark vessel ------------------------

Price_Bench = 371;                        %Price benchmark vessel (million EUR)
Cap_Bench = 11000;                        %Capacity of benchmark vessel (ton of H2)
Cap_LH2 = 1000:100:14000;                 %Range of capacities tested
Price_LH2_m = [];                         %Creating matrix

irr = 0.15;                               %Internal rate of return
OPEX_ship = Price_Bench*OPEX_rate_ship;   %OPEX is equal to 4% of CAPEX
charter_rate = 1:2:250;                   %Chater rates range to test
Life_ship = 5:1:25;                       %Startin at 5 to account for the design and construction period
NPV_charter = [];                         %Creating matrix

for c = charter_rate                      %Testing what charter rate gets a internal rate of return over the ships operation
  OPEX_Charter_irr = 0;                   %Restarting value for each iteraction
  for l = Life_ship
    OPEX_Charter_irr = OPEX_Charter_irr + (OPEX_ship-c)/((1+irr)^l);
  end
  NPV_charter = [NPV_charter,(Price_Bench+OPEX_Charter_irr)];
end

Point = find((abs(NPV_charter))<5);       %Checking what point is NPV closest to zero
if mod(length(Point),2)==1                %If captures more than one point, make sure to get the medium one
  Point = median(Point);
else                                      %If captures more than one point, make sure to get the medium one
  Point = Point((length(Point)/2)+1);
end

Ship_cap_freight_ratio = charter_rate(Point)/Price_Bench; %Relation between ship building cost and the charter rate
                                                          %according to a certain given internal rate of return (irr)
%---------------------------------------------------------------------------------------------------------------------

for i = Dist
    LCOH2_p_m = [];                       %Creating matrix

    LCOH2_sto_m = [];                     %Creating matrix
    LCOH2_ship_m = [];                    %Creating matrix
    LCOH2_liq = 0;

    for k = Power
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
        CAPEX_p_km = D*857000*0.9/25; %((4000000*(Diameter^2))+(598600*Diameter)+329000)/3 %Capital cost per km (EUR/km]
        CAPEX_com = (n_pump*2000000)*0.9;
        CAPEX_p = ((i*CAPEX_p_km) + CAPEX_com); %*((1+dr)^7);    %CAPEX multiplied by distance brought to present value 2022 (source 2010)

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


        %---Calculating the NPV of OPEX for pipeline and ship costs -------------------------------------%

        OPEX_p = 0;            %Initiating and restarting value for calculating NPV for next distance
        OPEX_liq = 0;          %Initiating and restarting value for calculating NPV for next distance
        OPEX_ship = 0;         %Initiating and restarting value for calculating NPV for next distance
               %Initiating and restarting value for calculating NPV for next distance
        for j = 1:1:40
          OPEX_p = OPEX_p + (CAPEX_p*OPEX_rate_p + (n_pump*408000*0.9))/((1+dr)^j);               %NPV of OPEX for pipeline
        endfor

        for j = L
            OPEX_ship = OPEX_ship + (CAPEX_stor*OPEX_rate_stor+...            %OPEX storage
                                     CAPEX_ship*OPEX_rate_ship+...            %OPEX tanker
                                     Fuel_cost)/((1+dr)^j);
        end

        LCOH2_p = (CAPEX_p+OPEX_p)/Levelized_prod;           %Levelized cost of hydrogen for pipeline transport[EUR/kgH2]
        LCOH2_p_m = [LCOH2_p_m,LCOH2_p];

        LCOH2_ship = (CAPEX_stor+OPEX_ship)/Levelized_prod;  %Levelized cost of hydrogen for ship transport (storage + ship + freight) [EUR/kgH2]
        LCOH2_ship_m = [LCOH2_ship_m,LCOH2_ship];  %Levelized cost of hydrogen for ship transport (Liquefaction + storage + ship + freight) [EUR/kgH2]


        if min([LCOH2_p,LCOH2_ship]) == LCOH2_p
          ZZ(find(Power==k),find(Dist==i))=1;
        elseif min([LCOH2_p,LCOH2_ship]) == LCOH2_ship
          ZZ(find(Power==k),find(Dist==i))=2;
        endif
    end
end

%subplot (2,2,4)
pcolor (Dist, Power, ZZ)
colormap ([0.4660, 0.6740, 0.1880; 0, 0.4470, 0.7410]);
set (gca,'fontsize',14,'fontname','times')
xlabel('Offshore distance (km)','fontname','times')
ylabel('Electrolyzer capacity (MW)','fontname','times')
%legend('Pipeline','HVDC','LH_tanker')
%legend ('Pipeline','fontsize',[18],'Location','Northeast')
%legend ('HVDC','fontsize',[18],'Location','Northeast')
%legend ('LH_{2} tanker','fontsize',[18],'Location','Northeast')

t = cputime

