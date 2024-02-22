clear all
%Inputs
Power = 500:50:3500;                     %Electrolyzer plant size [MW]
Dist = 100:50:3000;              %Distance [km]
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
OPEX_rate_HVDC = 0.01;                    %OPEX = 1% of CAPEX (Van Eeckhout et al 2009)
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
    LCOH2_liq_m = [];                     %Creating matrix
    LCOH2_liq_m = [];                     %Creating matrix
    LCOH2_sto_m = [];                     %Creating matrix
    LCOH2_ship_m = [];                    %Creating matrix
    LCOH2_HVDC_m = [];                    %Creating matrix
    Power_HVDC = [];                      %Creating matrix
    LCOH2_liq = 0;

    for k = Power
        %Control variables
        Prod_h = k*Elec_prod*CF*3600;      %Amount of hydrogen produced per hour [kg/h]
        Prod_d = Prod_h*24;                %Amount of hydrogen produced per day [kg/d]
        Prod_y = Prod_d*350;               %Amount of hydrogen produced per year [kg/y]
        Prod_l = Prod_y*max(L);            %Amount of hydrogen produced over the plant lifetime [kg]
        Levelized_prod = 0;                %Initiating variable
        for j = L
            Levelized_prod = Levelized_prod + (Prod_y/((1+dr)^j));  %NPV of production
        end

%======================================= PIPELINES =======================================%

        Q =k*Elec_prod*CF/Urate_p;             %Peak hydrogen mass flow rate [kg/s]
        H_den = 8;                             %Mass density of hydrogen [kg/m3]
        vel = 15;                              %Average fluid velocity [m/s]

        Diameter = sqrt((4*Q)/(H_den*vel*pi)); %Diameter in [m]
        if Diameter<0.1
        Diameter=0.1;                          %Minimum diammeter
        end

        CAPEX_p_km = ((4000000*(Diameter^2))+(598600*Diameter)+329000); %Capital cost per km (EUR/km]
        CAPEX_p = (i*CAPEX_p_km)*((1+dr)^12)*2;   %CAPEX multiplied by distance brought to present value 2022 (source 2010)

%======================================== SHIP ========================================%

        %--------------------  CAPEX liquefaction  --------------------%
        Cap_l = Prod_d/1000;                  %Design capacity [ton/day]
        Liq_max = 200;                        %Liquifier maximum capacity per unity (200,000 kg/day) [ton/day]
        N = Cap_l/Liq_max;                    %Number of liquefier (one liquefier can supply up to 200 ton/day)
        I = 1.16*0.9;                         %Cost index for cost estimates of 2016 source: DOE (2019)

        if N>1  %In case the need of more than one liquifier
            %Installed liquifier CAPEX brought to present value (source 2016)
            N = floor(N);
            CAPEX_liq = (((Liq_max^0.8*N*5.6*I)+((Cap_l-(Liq_max*N))^0.8*5.6*I))*1000000)*((1+dr)^6);

        else    %Just one liquified required
            %Installed liquifier CAPEX brought to present value 2022 (source 2016)
            CAPEX_liq = (Cap_l^0.8*5.6*I*1000000)*((1+dr)^6);
        end

        %-------------------- CAPEX Storage  --------------------%
        Ship_cap = Prod_h*2*i/V_ship;
        Storage_cap = Ship_cap*2;

        CAPEX_stor = Storage_cap*Storage_cost*((1+dr)^3);

        %--------------------  Charter rates  --------------------%

        CAPEX_ship = (Price_Bench*(Ship_cap/(Cap_Bench*1000))^0.623)*1000000;
        Charter_price_LH2 = CAPEX_ship*Ship_cap_freight_ratio;

%======================================== HVDC cable ========================================%
        if k<=2500
          CAPEX_HVDC = (444000*exp(0.00025*k/0.6))*i*((1+dr)^7)...  %Equation describing cost of HVDC offshore cables (Appendix E - National Grid Eso 2015)
                        + 600000*i*((1+dr)^7)...                    %Installation costs of HVDC offshore cables (Appendix E - National Grid Eso 2015)
                        + (334000*k/0.6-67120000)*((1+dr)^5);   %Converter costs of VSC HVDC interconnectors (Hartel et al 2017)
          Power_HVDC = [Power_HVDC, k];
        endif

        %---Calculating the NPV of OPEX for pipeline and ship costs -------------------------------------%

        OPEX_p = 0;            %Initiating and restarting value for calculating NPV for next distance
        OPEX_liq = 0;          %Initiating and restarting value for calculating NPV for next distance
        OPEX_ship = 0;         %Initiating and restarting value for calculating NPV for next distance
        OPEX_HVDC = 0;         %Initiating and restarting value for calculating NPV for next distance

        for j = L
            OPEX_p = OPEX_p + (CAPEX_p*OPEX_rate_p)/((1+dr)^j);               %NPV of OPEX for pipeline
            OPEX_liq = OPEX_liq + (CAPEX_liq*OPEX_rate_liq)/((1+dr)^j);       %NPV of OPEX for liquefaction
            OPEX_ship = OPEX_ship + (CAPEX_stor*OPEX_rate_stor+...            %OPEX storage
                                     CAPEX_ship*OPEX_rate_ship+...            %OPEX tanker
                                     Charter_price_LH2)/((1+dr)^j);           %Charter rates
            OPEX_HVDC = OPEX_HVDC + (CAPEX_HVDC*OPEX_rate_HVDC)/((1+dr)^j);   %OPEX HVDC
        end

        LCOH2_p = (CAPEX_p+OPEX_p)/Levelized_prod;           %Levelized cost of hydrogen for pipeline transport[EUR/kgH2]
        LCOH2_p_m = [LCOH2_p_m,LCOH2_p];

        LCOH2_liq = max(LCOH2_liq,(CAPEX_liq+OPEX_liq)/Levelized_prod);  %Levelized cost of hydrogen related to liquefaction [EUR/kgH2]
        LCOH2_liq_m = [LCOH2_liq_m,LCOH2_liq];

        LCOH2_ship = (CAPEX_stor+OPEX_ship)/Levelized_prod;  %Levelized cost of hydrogen for ship transport (storage + ship + freight) [EUR/kgH2]
        LCOH2_ship_m = [LCOH2_ship_m,LCOH2_ship+LCOH2_liq];  %Levelized cost of hydrogen for ship transport (Liquefaction + storage + ship + freight) [EUR/kgH2]

        if k<=2500
          LCOH2_HVDC = (CAPEX_HVDC+OPEX_HVDC)/Levelized_prod;
          LCOH2_HVDC_m = [LCOH2_HVDC_m, LCOH2_HVDC];
        endif

        if min([LCOH2_p,LCOH2_ship+LCOH2_liq,LCOH2_HVDC]) == LCOH2_p
          ZZ(find(Power==k),find(Dist==i))=1;
        elseif min([LCOH2_p,LCOH2_ship+LCOH2_liq,LCOH2_HVDC]) == LCOH2_HVDC
          ZZ(find(Power==k),find(Dist==i))=2;
        elseif min([LCOH2_p,(LCOH2_ship+LCOH2_liq),LCOH2_HVDC]) == LCOH2_ship+LCOH2_liq
          ZZ(find(Power==k),find(Dist==i))=3;
        endif
    end
end

%subplot (2,2,4)
pcolor (Dist, Power, ZZ)
colormap ([0.4660, 0.6740, 0.1880; 0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410]);
set (gca,'fontsize',14,'fontname','times')
xlabel('Offshore distance (km)','fontname','times')
ylabel('Electrolyzer capacity (MW)','fontname','times')
%legend ('Pipeline','fontsize',[18],'Location','Northeast')
%legend ('LH_{2} tanker','fontsize',[18],'Location','Northeast')
%legend ('HVDC','fontsize',[18],'Location','Northeast')

t = cputime

