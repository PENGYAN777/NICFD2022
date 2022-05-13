% shock loss calculation

%% PREAMBLE
close all
clear all

%% DATA
% PT    = 11.12;                                                              % [bar]
% TT    = 230.3;                                                              % [°C]
% Pfs   = 3.6;                                                                % [bar]
% fluid = 'MM';
% 
PT    = 8.3894;                                                                % [bar]
TT    = 18.4463;                                                                 % [°C]
Pfs   = 1.6765;                                                                % [bar]
fluid = 'nitrogen';


%% CoolProp
MM    =  PropsSI('molar_mass',fluid);
R     = 8.314472/MM;                                                          % costante universale dei gas [J/kg K]

%% IDEAL GAS SHOCK RELATIONS with free stream conditions
[DPt_IG,gamma_IG,rho2rho1_IG] = fun_IGshock(PT,TT,Pfs,fluid);
Pt2_IG                        = PT - DPt_IG;
Y_IG                          = (PT - Pt2_IG)./(PT - Pfs)*100;

%% non IG shock relation
[~,~,~,~,Pt2_nIG,Tt2_nIG,~] = fun_shockIt_rho2(PT,TT,Pfs,rho2rho1_IG,fluid);
DPt_nIG_Pfs                 = PT - Pt2_nIG;
Y_nIG                       = (PT - Pt2_nIG)./(PT - Pfs)*100;


%%
% ------------------------------------------------------------------------
% ---------------             THE END                ---------------------
% ------------------------------------------------------------------------


%% FUNCTIONS

% SHOCK nonIG, iterative in terms of rho2: inputs are PT[bar], TT[°C], Pfs[bar]=pre-shock static pressure (CoolProp)
function [P2,RHO2,U2,H2,PT2,TT2,violated] = fun_shockIt_rho2(PT,TT,Pfs,rho2rho1_IG,fluid)

    P1   = Pfs;
    S1   = PropsSI('S','P',PT*10^5,'T',TT+273.15,fluid)/10^3;               % [kJ/kg/K]
    HT1  = PropsSI('H','P',PT*10^5,'T',TT+273.15,fluid)/10^3;               % [kJ/kg]
    H1   = PropsSI('H','P',P1*10^5,'S',S1*10^3,fluid)/10^3;                 % [kJ/kg]
    U1   = sqrt(2*(HT1-H1)*10^3);
    RHO1 = PropsSI('D','P',P1*10^5,'S',S1*10^3,fluid);
       
    function err = Rank_hug(x)
        U2  = RHO1.*U1./x;
        P2  = (P1*10^5 + RHO1.*U1.^2 - x.*U2.^2)/10^5;                      % [bar]
        H2  = (H1*10^3 + 1/2*U1.^2  - 1/2*U2.^2)/10^3;                      % [kJ/kg]
        err = (H2 - PropsSI('H','P',P2*10^5,'D',x,fluid)/10^3)*10^8;
    end

    options                      = optimoptions('fsolve','Display','final');
    options.FunctionTolerance    = 1e-8;
    options.OptimalityTolerance  = 1e-8; 
    options.StepTolerance        = 1e-8; 
    options.FiniteDifferenceType = 'central';
    
    x0    = RHO1.*rho2rho1_IG; % uso IG per inizializzare
    [x,~] = fsolve(@Rank_hug,x0,options);
    RHO2  = x;

    S2  = PropsSI('S','P',P2*10^5,'D',RHO2,fluid)/10^3;                     % [kJ/kg/K]
    HT2 = (H2*10^3 + U2.^2/2)/10^3;
    PT2 = PropsSI('P','H',HT2*10^3,'S',S2*10^3,fluid)/10^5;                 % [bar]
    TT2 = PropsSI('T','H',HT2*10^3,'S',S2*10^3,fluid)-273.15;               % [°C]

    % check bal
    cons_mass = RHO1.*U1 - RHO2.*U2;
    cons_mom  = P1*10^5 + RHO1.*U1.^2 - P2*10^5 - RHO2.*U2.^2;
    cons_en   = H1*10^3 + 1/2*U1.^2 - H2*10^3 - 1/2*U2.^2;
    
   violated = 0; % []
   if any(cons_mass > 10^-8)
       disp('Watch out! Mass conservation across shock possibly violated')
%        violated = find(cons_mass > 10^-8);
         violated = 1;
   end
   if any(cons_mom  > 10^-7)
       disp('Watch out! Momentum conservation across shock possibly violated')
   end
   if any(cons_en   > 10^-6)
       disp('Watch out! Energy conservation across shock possibly violated')
   end
    
end

% SHOCK IG: inputs are PT[bar], TT[°C], Pfs[bar]=pre-shock static pressure (CoolProp)
function [DPT_IG,gamma_IG,R2R1_IG] = fun_IGshock(PT,TT,Pfs,fluid)
    
    p_idealgas = 10e-7*ones(size(PT,1),1);

    P1       = Pfs;
    HT1      = PropsSI('H','P',PT*10^5,'T',TT+273.15,fluid)/10^3;           % [kJ/kg]
    S1       = PropsSI('S','P',PT*10^5,'T',TT+273.15,fluid)/10^3;           % [kJ/kg/K]
    H1       = PropsSI('H','P',P1*10^5,'S',S1*10^3,fluid)/10^3;             % [kJ/kg]
    U1       = sqrt(2*(HT1-H1)*10^3);
    A1       = PropsSI('A','P',P1*10^5,'S',S1*10^3,fluid) ;                 % [kJ/kg]
    M1       = U1./A1;
    T1       = PropsSI('T','P',P1*10^5,'S',S1*10^3,fluid)-273.15;           % [°C]
    cp_IG    = PropsSI('C','P',p_idealgas*10^5,'T',T1+273.15,fluid)/10^3;   % [kJ/kg/K]
    cv_IG    = PropsSI('O','P',p_idealgas*10^5,'T',T1+273.15,fluid)/10^3;   % [kJ/kg/K]
    gamma_IG = cp_IG./cv_IG;

    PT2_PT1  = (((gamma_IG+1).*M1.^2)./((gamma_IG-1).*M1.^2+2)).^(gamma_IG./(gamma_IG-1))...
                     .*((gamma_IG+1)./((2.*gamma_IG.*M1.^2-(gamma_IG-1)))).^(1./(gamma_IG-1));
    PT2_IG   = PT2_PT1.*PT;
    DPT_IG   = PT - PT2_IG;

    R2R1_IG  = ((gamma_IG+1).*M1.^2)./((gamma_IG-1).*M1.^2+2);              % per inizializzare proc. it. 
    
end

