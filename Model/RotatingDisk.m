%% Matlab Initilizations

clear;clc;close all;
lc = 0:0.001:1;
res_lc = 0:0.01:1;

%% Defining the material parameters

syms rho sig_o nu E eta

%Values are taken for Structural Steel(Non-Linear) in SI Units
nu_val = 0.3;
eta_val = 5.8785;

%Values used in the reference paper
%nu_val = 1/3; 
%eta_val = 400; 
%eta_val = 8e-4; 

sig_o_val = 2.5e8;
rho_val = 7850;
E_val = 2e11;

H = (eta*sig_o)/E;
W = sqrt(H/(1+H));
A = 1-nu+(1/(2*H));

%% Defining the Disk Parameters 

syms b r r1 r2 x x1 x2 bet omg Omg 

b_val = 0.5;

%omg_val = 591.3980;      % Value used in the reference paper for x2 = 0.5
%omg_val = 742.5364;      % Value used in the reference paper for x2 = 1
%omg_val = 1.0707*1.0e+3; % Value used in the reference paper for x2 = and H = 10e-6 (Incorrect Value)
%omg_val = 618;
%omg_val = 500;           % Value for elastic
omg_val = 600;    

disp('Input angular velocity is (in rad/s):')
disp(omg_val)

Omg_eqn =  Omg == sqrt((rho*(omg^2)*(b^2))/sig_o);              

%% Elastic Region

syms u_e sig_r_e sig_t_e

ue_eqn = u_e == ((1 - nu)/8)*(Omg^2)*x*(3+nu-((1+nu)*x^2));
sige_eqn = sig_r_e == ((3+nu)/8)*(Omg^2)*(1-(x^2));
sigte_eqn = sig_t_e == (1/8)*(Omg^2)*(3+nu-((1+(3*nu))*x^2));

%double(solve(subs(Omg_eqn,[Omg,rho,b,sig_o],[3.00001,rho_val,b_val,sig_o_val]),omg))

%% Plastic Region 

%Terms for equation 3.25
T1 = ((1+(3*nu))/(9+(8*H)))*(bet^3);
T2 = H/(9+(8*H));     
T3 = 1-(2*(1+(3*nu))*T2);
T4 = (1/W) - (W*A);
T5 = ((1-(3*nu))/9)-A;
T6 = T5*T2*(1+(3*nu));  
T7 = (1/9)+T6+(A/4);

bet_eqn = T1+((bet^W)*(((-1/3)*T3*(T4+nu))-((1+W)*T7)))+((bet^(-W))*(((1/3)*T3*(T4-nu))-((1-W)*T7))) == 0;            % Equation 3.25 from the paper
A_ = subs(bet_eqn,[eta,sig_o,E,nu],[eta_val,sig_o_val,E_val,nu_val]);                                                 % Substitute the values of all parameters to form a equation in beta
bet_val = double(vpasolve(A_,bet,1000));                                                                              % Solve the equation for value of beta 

% Just for reference
%{
disp('The value of beta is:')
disp(bet_val)
%}

%Terms for equation 3.23
P1 = ((3+nu)/4)*(1/(1+(x2^2)));
P2 = ((1+(3*nu))/8)*((1-(x2^2))/(1+(x2^2)));
P3 = (3+nu)/8;
P4 = 1+(3*nu);
P5 = (1/W)-(1/3); 
P6 = (1/W)+(1/3);
P7 = P4*T2*P5;
P8 = P4*T2*P6;
B1 = ((-P1)+((P2+P3-(1/3)+P7)*(x2^2)))*(bet^(3+W))*(-(1+W));
B2 = (P1+((-P2-P3+(1/3)+P8)*(x2^2)))*(bet^(3-W))*(1-W);
B3 = (2/3)*T3*(x2^2);
B4 = (1+W)*(bet^(3+W));
B5 = (1-W)*(bet^(3-W));
B6 = (2/(1+(x2^2)))*(B4+B5);

x2_eqn = (Omg^2)*(B1+B2-B3) == B6;                                                                                    % Equation 3.23 from the paper
B = subs(x2_eqn,[eta, sig_o, E, nu, bet], [eta_val, sig_o_val, E_val, nu_val, bet_val]);                              % Substituting the value of all parameters to form a equation in x2
Omg_val_eq = solve(B,Omg);                                                                                            % Solve the equation in terms of x2 for Omg
figure(1)
fplot(Omg_val_eq(2),[0 1])                                                                                            % Plot Omg vs x2. Compare with Fig.1 in the paper
title('Variation of Non-Dimensional Angular Velocity with x')
xlabel('x');
ylabel('Omg');
ylim([1 2.5]);
                                                                                                                                                      
Omg_val_arr = double(subs(solve(B,Omg),x2,lc));                                                                       % Get the values of Omg_val                                    
Omg_val_arr = Omg_val_arr(2,:);                                                                                       % Discarding the negative values

% Just for reference
%{
omg_2_val = subs(solve(Omg_eqn, omg),[Omg,rho,b,sig_o],[Omg_val(2),rho_val,b_val,sig_o_val]);                         % Solve for the value of omg2 in terms of x2
figure(2)
fplot(omg_2_val(1),[0 1])                                                                                             % Plot omg vs x2
title('Variation of Angular Velocity with x')
xlabel('x');
ylabel('omg');
%}

C = subs(Omg_eqn,[rho,b,sig_o],[rho_val,b_val,sig_o_val]);                                                            % Substituting the value of all parameters except and forming an equation in Omg
omg_val_arr = double(subs(solve(C,omg),Omg,Omg_val_arr));                                                             % Solving for omg for evry value of Omg obtained earlier
omg_val_arr = omg_val_arr(2,:);                                                                                       % Discarding the negative values
disp('Angular velocity (in rad/s) at which yielding starts at the centre is')
disp(omg_val_arr(1))
disp('Angular velocity (in rad/s) at which the disc becomes fully plastic is')
disp(omg_val_arr(end))

%% Calculating Stresses, Displacement and Strains

Omg_val = double(subs(solve(Omg_eqn, Omg),[rho, omg, b, sig_o],[rho_val, omg_val, b_val, sig_o_val]));                % Calculating the value of Omg for further calculations

% Checking the value of angular velocity to determine which regions are formed
if omg_val < omg_val_arr(1)                                                                                           % Checking if the disc is fully elastic                                                                                                                                 
    
    u_e_arr = double(subs(solve(subs(ue_eqn,[Omg,nu],[Omg_val,nu_val]),u_e),x,res_lc));                               % Making an array of displacements
    [u_e_min, ind_u_e_min] = min(u_e_arr);                                                                            % Getting minimum displacement
    [u_e_max, ind_u_e_max] = max(u_e_arr);                                                                            % Getting maximum displacement
    fprintf('%d (in mm) is the maximum displacement at %d (in mm)\n',1000*u_e_max*((b_val*sig_o_val)/E_val),res_lc(ind_u_e_max)*b_val*1000)                       
    fprintf('%d (in mm) is the minimum displacement at %d (in mm)\n\n',1000*u_e_min*((b_val*sig_o_val)/E_val),res_lc(ind_u_e_min)*b_val*1000)
    
    sig_r_e_arr = double(subs(solve(subs(sige_eqn,[Omg,nu],[Omg_val,nu_val]),sig_r_e),x,res_lc));                     % Making an array of radial stresses
    [sig_r_e_min, ind_r_e_min] = min(sig_r_e_arr);                                                                    % Getting minimum radial stresses
    [sig_r_e_max, ind_r_e_max] = max(sig_r_e_arr);                                                                    % Getting maximum radial stresses
    fprintf('%d (in Pa) is the maximum radial stress at %d (in mm)\n',sig_r_e_max*sig_o_val,res_lc(ind_r_e_max)*b_val*1000)
    fprintf('%d (in Pa) is the minimum radial stress at %d (in mm)\n\n',sig_r_e_min*sig_o_val,res_lc(ind_r_e_min)*b_val*1000)
    
    sig_t_e_arr = double(subs(solve(subs(sigte_eqn,[Omg,nu],[Omg_val,nu_val]),sig_t_e),x,res_lc));                    % Making an array of tangential stresses
    [sig_t_e_min, ind_t_e_min] = min(sig_t_e_arr);                                                                    % Getting minimum tangential stresses
    [sig_t_e_max, ind_t_e_max] = max(sig_t_e_arr);                                                                    % Getting maximum tangential stresses
    fprintf('%d (in Pa) is the maximum tangential stress at %d (in mm)\n',sig_t_e_max*sig_o_val,res_lc(ind_t_e_max)*b_val*1000)
    fprintf('%d (in Pa) is the minimum tangential stress at %d (in mm)\n\n',sig_t_e_min*sig_o_val,res_lc(ind_t_e_min)*b_val*1000)
    
    figure(3)
    fplot(solve(subs(ue_eqn,[Omg,nu],[Omg_val,nu_val]),u_e),[0 1])                                                    % Solving the displacement equation and plotting the displacement variation with radius
    hold on
    fplot(solve(subs(sige_eqn,[Omg,nu],[Omg_val,nu_val]),sig_r_e),[0 1])                                              % Solving the radial stress equation and plotting the radial stress variation with radius
    hold on
    fplot(solve(subs(sigte_eqn,[Omg,nu],[Omg_val,nu_val]),sig_t_e),[0 1])                                             % Solving the tangential stress equation and plotting the tangential stress variation with radius
    title('Variation of Radial & Tangential Stress as well as Displacement with radius when the disk is completely Elastic')
    xlabel('x');
    legend('Displacement','Radial Stress','Tangential Stress');
else 
%elseif omg_val_arr(1) < omg_val < omg_val_arr(end)
    
    [Omg_val_err, ind] = min(abs(Omg_val_arr-Omg_val));
    %x2_val = 0.5; 
    x2_val = lc(ind);
    Omg_val = Omg_val_arr(ind);
    x1_val = round(double(x2_val/bet_val),2);
    fprintf('Plastic Region II starts at %d (in mm)\n',x1_val*b_val*1000);
    fprintf('Elastic Region starts at %d (in mm)\n',x2_val*b_val*1000);
    
    %Terms for equation 3.18, 3.19, 3.20 and 3.27
    P9 = 1/(1+(x2^2));
    P10 = P3*(Omg^2);
    
    syms C1 C3 C4 C5
    
    C5 = ((x2^2)*P9)*(1-((1/8)*(Omg^2)*(3+nu-(P4*(x2^2)))));                                                          % Equation 3.18 from the paper
    C3 = (x2^(1+W))*(((-P9)*(P10-1))+((1/2)*(P2+P3-(1/3)+P7)*(Omg^2)*(x2^2)));                                        % Equation 3.19 from the paper
    C4 = (x2^(1-W))*((P9*(P10-1))+((1/2)*(-P2-P3+(1/3)+P8)*(Omg^2)*(x2^2)));                                          % Equation 3.20 from the paper
    C1 = ((-C3)*(x1^(-(1+W))))+(C4*(x1^(-(1-W))))+((1/3)*((1/2)-(P4*T2))*(Omg^2)*(x1^2))+1;                           % Equation 3.27 from the paper 
    
    %Terms for equation 3.30, 3.31 and 3.32
    P11 = (1/W) + nu;
    P12 = (1/W) - nu;
    P13 = 1 - (9*(nu^2));
    P14 = P13*T2;
    
    %Terms for equation 3.33 and 3.35
    P15 = (1-nu)/8;
    P16 = (1+nu)*(x^2);
    
    % Plastic Region I
    syms u_p1 sig_r_p1 sig_t_p1 str_r_p1 str_t_p1 str_z_p1
    
    up1_eqn = u_p1 == ((1-nu+(1/(2*H)))*x*(C1-((1/4)*(Omg^2)*(x^2)))) - (x*(1/(2*H)));                                % Equation 3.28 from the paper
    sigp1_eqn = sig_r_p1 == C1 - ((1/2)*(Omg^2)*(x^2));                                                               % Equation 3.29 from the paper
    sigtp1_eqn = sig_t_p1 == C1 - ((1/2)*(Omg^2)*(x^2));                                                              % Equation 3.29 from the paper
    strp1_eqn = str_r_p1 == ((1/(2*H))*(C1-1))-((1/4)*(1-nu+(3/(2*H)))*(x^2)*(Omg^2));                                % Equation 3.36 from the paper
    strtp1_eqn = str_t_p1 == ((1/(2*H))*(C1-1))+((1/4)*(1-nu-(1/(2*H)))*(x^2)*(Omg^2));                               % Equation 3.37 from the paper
    strzp1_eqn = str_z_p1 == (-(1/H)*(C1-1)) + ((1/(2*H))*(x^2)*(Omg^2));                                             % Equation 3.38 from the paper
    
    u_p1_arr = double(subs(solve(subs(up1_eqn,[nu eta sig_o E x1 x2 Omg],[nu_val eta_val sig_o_val E_val x1_val x2_val Omg_val]),u_p1),x,0:0.001:x1_val));                % Making an array of displacements in Plastic Region I                 
    sig_r_p1_arr = double(subs(solve(subs(sigp1_eqn,[nu eta sig_o E x1 x2 Omg],[nu_val eta_val sig_o_val E_val x1_val x2_val Omg_val]),sig_r_p1),x,0:0.001:x1_val));      % Making an array of radial stresses in Plastic Region I
    sig_t_p1_arr = double(subs(solve(subs(sigtp1_eqn,[nu eta sig_o E x1 x2 Omg],[nu_val eta_val sig_o_val E_val x1_val x2_val Omg_val]),sig_t_p1),x,0:0.001:x1_val));     % Making an array of tangential stresses in Plastic Region I
    str_r_p1_arr = double(subs(solve(subs(strp1_eqn,[nu eta sig_o E x1 x2 Omg],[nu_val eta_val sig_o_val E_val x1_val x2_val Omg_val]),str_r_p1),x,0:0.001:x1_val));      % Making an array of radial strains in Plastic Region I
    str_t_p1_arr = double(subs(solve(subs(strtp1_eqn,[nu eta sig_o E x1 x2 Omg],[nu_val eta_val sig_o_val E_val x1_val x2_val Omg_val]),str_t_p1),x,0:0.001:x1_val));     % Making an array of tangential strains in Plastic Region I
    str_z_p1_arr = double(subs(solve(subs(strzp1_eqn,[nu eta sig_o E x1 x2 Omg],[nu_val eta_val sig_o_val E_val x1_val x2_val Omg_val]),str_z_p1),x,0:0.001:x1_val));     % Making an array of planar strains in Plastic Region I
    
    % Plastic Region II
    syms u_p2 sig_r_p2 sig_t_p2 str_r_p2 str_t_p2 str_z_p2
    
    up2_eqn = u_p2 == (P11*C3*(x^(-W)))+(P12*C4*(x^W))-((1/9)*(1+P14)*(Omg^2)*(x^3))+((1-nu)*x);                                             % Equation 3.30 from the paper
    sigp2_eqn = sig_r_p2 == ((-C3)*(x^(-(1+W))))+(C4*(x^(-(1-W))))-((1/3)*(1+(P4*T2))*(Omg^2)*(x^2))+1;                                      % Equation 3.31 from the paper
    sigtp2_eqn = sig_t_p2 == (W*C3*(x^(-(1+W))))+(W*C4*(x^(-(1-W))))-(P4*T2*(Omg^2)*(x^2))+1;                                                % Equation 3.32 from the paper
    strp2_eqn = str_r_p2 == (1/E)*(solve(sigtp2_eqn,sig_t_p2)-(nu*(solve(sigp2_eqn,sig_r_p2))));                                             % Equation 2.3 from the paper
    strtp2_eqn = str_t_p2 == (((1-(W^2))/W)*((C3*(x^(-(1+W))))+(C4*(x^(-(1-W))))))-(((1+(3*nu))/9)*(1-((8*H)/(9+(8*H))))*(Omg^2)*(x^2));     % Equation 3.39 from the paper
    strzp2_eqn = str_z_p2 == -((((1-(W^2))/W)*((C3*(x^(-(1+W))))+(C4*(x^(-(1-W))))))-(((1+(3*nu))/9)*(1-((8*H)/(9+(8*H))))*(Omg^2)*(x^2)));  % Equation 3.39 from the paper
    
    u_p2_arr = double(subs(solve(subs(up2_eqn,[nu,sig_o,E,eta,Omg,x2],[nu_val,sig_o_val,E_val,eta_val,Omg_val,x2_val]),u_p2),x,x1_val:0.001:x2_val));                     % Making an array of displacements in Plastic Region II
    sig_r_p2_arr = double(subs(solve(subs(sigp2_eqn,[nu,sig_o,E,eta,Omg,x2],[nu_val,sig_o_val,E_val,eta_val,Omg_val,x2_val]),sig_r_p2),x,x1_val:0.001:x2_val));           % Making an array of radial stresses in Plastic Region II
    sig_t_p2_arr = double(subs(solve(subs(sigtp2_eqn,[nu,sig_o,E,eta,Omg,x2],[nu_val,sig_o_val,E_val,eta_val,Omg_val,x2_val]),sig_t_p2),x,x1_val:0.001:x2_val));          % Making an array of tangential stresses in Plastic Region II
    str_r_p2_arr = double(subs(solve(subs(strp2_eqn,[nu eta sig_o E x1 x2 Omg],[nu_val eta_val sig_o_val E_val x1_val x2_val Omg_val]),str_r_p2),x,x1_val:0.001:x2_val));      % Making an array of radial strains in Plastic Region II
    str_t_p2_arr = double(subs(solve(subs(strtp2_eqn,[nu eta sig_o E x1 x2 Omg],[nu_val eta_val sig_o_val E_val x1_val x2_val Omg_val]),str_t_p2),x,x1_val:0.001:x2_val));     % Making an array of tangential strains in Plastic Region II
    str_z_p2_arr = double(subs(solve(subs(strzp2_eqn,[nu eta sig_o E x1 x2 Omg],[nu_val eta_val sig_o_val E_val x1_val x2_val Omg_val]),str_z_p2),x,x1_val:0.001:x2_val));     % Making an array of planar strains in Plastic Region II
   
    % Elastic Region
    syms up_e sigp_r_e sigp_t_e
    
    up_e_eqn = up_e == (C5*(((1+nu)*(1/x))+((1-nu)*x)))+(P15*(Omg^2)*x*(3+nu-P16));                                   % Equation 3.33 from the paper
    sigp_e_eqn = sigp_r_e == ((-C5)*((1/(x^2))-1))+(P3*(Omg^2)*(1-(x^2)));                                            % Equation 3.34 from the paper
    sigtp_e_eqn = sigp_t_e == (C5*((1/(x^2))+1))+((1/8)*(Omg^2)*(3+nu-(P4*(x^2))));                                   % Equation 3.35 from the paper
    
    up_e_arr = double(subs(solve(subs(up_e_eqn,[nu,sig_o,E,eta,Omg,x2],[nu_val,sig_o_val,E_val,eta_val,Omg_val,x2_val]),up_e),x,x2_val:0.001:1));                         % Making an array of displacements in Elastic Region
    sigp_r_e_arr = double(subs(solve(subs(sigp_e_eqn,[nu,sig_o,E,eta,Omg,x2],[nu_val,sig_o_val,E_val,eta_val,Omg_val,x2_val]),sigp_r_e),x,x2_val:0.001:1));               % Making an array of radial stresses in Elastic Region
    sigp_t_e_arr = double(subs(solve(subs(sigtp_e_eqn,[nu,sig_o,E,eta,Omg,x2],[nu_val,sig_o_val,E_val,eta_val,Omg_val,x2_val]),sigp_t_e),x,x2_val:0.001:1));              % Making an array of tangential stresses in Elastic Region 
    
    % Combining all the values
    up = [u_p1_arr(1:end-1) (u_p1_arr(end)+u_p2_arr(1))/2 u_p2_arr(2:end-1) (u_p2_arr(end)+up_e_arr(1))/2  up_e_arr(2:end)];                                              % Forming a combined array of displacements
    sigp_r = [sig_r_p1_arr(1:end-1) (sig_r_p1_arr(end)+sig_r_p2_arr(1))/2 sig_r_p2_arr(2:end-1) (sig_r_p2_arr(end)+sigp_r_e_arr(1))/2 sigp_r_e_arr(2:end)];               % Forming a combined array of radial stresses
    sigp_t = [sig_t_p1_arr(1:end-1) (sig_t_p1_arr(end)+sig_t_p2_arr(1))/2 sig_t_p2_arr(2:end-1) (sig_t_p2_arr(end)+sigp_t_e_arr(1))/2 sigp_t_e_arr(2:end)];               % Forming a combined array of tangential stresses
    strp_r = [str_r_p1_arr(1:end-1) (str_r_p1_arr(end)+str_r_p2_arr(1))/2 str_r_p2_arr(2:end)];                                                                           % Forming a combined array of radial strains
    strp_t = [str_t_p1_arr(1:end-1) (str_t_p1_arr(end)+str_t_p2_arr(1))/2 str_t_p2_arr(2:end)];                                                                           % Forming a combined array of tangential strains
    strp_z = [str_z_p1_arr(1:end-1) (str_z_p1_arr(end)+str_z_p2_arr(1))/2 str_z_p2_arr(2:end)];                                                                           % Forming a combined array of planar strains
      
    X =  [0:0.001:x1_val (x1_val+0.001):0.001:x2_val (x2_val+0.001):0.001:1];                                                                                             % Forming an array for X-axis for stresses and displacements
    Xr = [0:0.001:x1_val (x1_val+0.001):0.001:x2_val];                                                                                                                    % Forming an array for X-axis for strains
    
    [up_min, ind_up_min] = min(up);
    [sigp_r_min, ind_sigp_r_min] = min(sigp_r);
    [sigp_t_min, ind_sigp_t_min] = min(sigp_t);
    [strp_r_min, ind_strp_r_min] = min(strp_r);
    [strp_t_min, ind_strp_t_min] = min(strp_t);
    [strp_z_min, ind_strp_z_min] = min(strp_z);
    
    [up_max, ind_up_max] = max(up);
    [sigp_r_max, ind_sigp_r_max] = max(sigp_r);
    [sigp_t_max, ind_sigp_t_max] = max(sigp_t);
    [strp_r_max, ind_strp_r_max] = max(strp_r);
    [strp_t_max, ind_strp_t_max] = max(strp_t);
    [strp_z_max, ind_strp_z_max] = max(strp_z);
    
    fprintf('\n%d (in mm) is the maximum displacement at %d (in mm)\n',1000*up_max*((b_val*sig_o_val)/E_val),X(ind_up_max)*b_val*1000)                       
    fprintf('%d (in mm) is the minimum displacement at %d (in mm)\n\n',1000*up_min*((b_val*sig_o_val)/E_val),X(ind_up_min)*b_val*1000)
    fprintf('%d (in Pa) is the maximum radial stress at %d (in mm)\n',sigp_r_max*sig_o_val,X(ind_sigp_r_max)*b_val*1000)
    fprintf('%d (in Pa) is the minimum radial stress at %d (in mm)\n\n',sigp_r_min*sig_o_val,X(ind_sigp_r_min)*b_val*1000)
    fprintf('%d (in Pa) is the maximum tangential stress at %d (in mm)\n',sigp_t_max*sig_o_val,X(ind_sigp_t_max)*b_val*1000)
    fprintf('%d (in Pa) is the minimum tangential stress at %d (in mm)\n\n',sigp_t_min*sig_o_val,X(ind_sigp_t_min)*b_val*1000)
    disp('Not Considering the Elastic Region')
    fprintf('%d  is the maximum radial strain at %d (in mm)\n',strp_r_max*(sig_o_val/E_val),X(ind_strp_r_max)*b_val*1000)
    fprintf('%d  is the minimum radial strain at %d (in mm)\n\n',strp_r_min*(sig_o_val/E_val),X(ind_strp_r_min)*b_val*1000)
    fprintf('%d  is the maximum tangential strain at %d (in mm)\n',strp_t_max*(sig_o_val/E_val),X(ind_strp_t_max)*b_val*1000)
    fprintf('%d  is the minimum tangential strain at %d (in mm)\n\n',strp_t_min*(sig_o_val/E_val),X(ind_strp_t_min)*b_val*1000)
    fprintf('%d  is the maximum planar strain at %d (in mm)\n',strp_z_max*(sig_o_val/E_val),X(ind_strp_z_max)*b_val*1000)
    fprintf('%d  is the minimum planar strain at %d (in mm)\n\n',strp_z_min*(sig_o_val/E_val),X(ind_strp_z_min)*b_val*1000)
                                                                                                 
    % Plotting the Stresses and Displacements
    figure(4)
    plot(X,up)
    hold on
    plot(X,sigp_r,'--')
    hold on
    plot(X,sigp_t,'g-.')
    title('Variation of Radial & Tangential Stress and Displacement with radius for Elastic Plastic Disk')
    %title('Variation of Radial & Tangential Stress and Displacement with radius for Completely Plastic Disk')
    xlabel('x');
    %xlim([0 0.01]);
    legend('Displacement','Radial Stress','Tangential Stress');
    
    figure(5)
    plot(Xr,strp_r)
    hold on
    plot(Xr,strp_t,'--')
    hold on
    plot(Xr,strp_z,'g-.')
    hold on
    title('Variation of Radial, Tangential and Planar Strains with radius for Elastic Plastic Disk')
    %title('Variation of Radial, Tangential and Planar Strains with radius for Completely Plastic Disk')
    xlabel('x');
    %xlim([0 0.01]);
    legend('Radial Strain','Tangential Strain','Planar Strain');    
end    