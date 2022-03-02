%% MAE 158 Drag Calculator
% Shashwat Sparsh
% ID: 21697442
% |Skin Friction Drag Coeff|
% Setup

hp_aircraft = 1;
T = 400;
%P = 1;  % via Table A.2
R = 1716;
rho = 0.0008754;

speedSound = sqrt(1.4*R*T)
V = [230:10:880]            %765;
M = V/speedSound            % Mach Number

Sref = 1000;

%CRexp = [];
%K = [];
%Cf = [];
%f = [];

%% Characteristic Lengths
MACexpW = 0;
MACexph = 0;
MACexpv = 0;
Lc = [MACexpW, MACexph, MACexpv, 16.2, 16.8];

%% Ratios
Lf = 119;
Df = 11;
ratios = [0.12, 0.09, 0.09, 0.06, Lf/Df, 5]  % Thickness and Fineness
sigma = [0.2, 0.35, 0.8, 1] % Taper

%% Sexp
Sexpw = 0;  % Defined Later on
Sexp = [Sexpw, 261, 161]

%% Swet
Swet = [0, 0, 0, 117, 0, 455]
% Wing

bW = 93.2;                      % Span
tcW = 0.18;
sweepangleW = 28;               % Sweep Angle
sigmaW = 0.2;                   % Taper Ratio: cT/cR
CRW = 17.8;                     % Root Chord
Coverage_wing = .17;            % Percent Covered
Rfuse = 11/2;

% Getting Swet %
SexpW = (1-Coverage_wing)*Sref;
SwetW = SWET(SexpW);

% Getting Skin Friction Coefficient %
CTW = sigmaW * CRW;
CRexp_wing = CREXP(CRW, CTW, Rfuse, bW);
MACexpW = MAC(CRexp_wing,CTW)
 
RNw = ReynoldsNumber(V, MACexpW);
Cf_w = CF(RNw);

% Getting Form Factor %
Kwing = Kairfoil(tcW, M, sweepangleW)

% Calculating f and adding to array %
fwing = F(Kwing, Cf_w, SwetW)
% Horizontal Tail

SexpH = 261;
tcH = 0.09;
sweepangleH = 31.6;             % Sweep Angle
sigmaH = 0.35;                  % Wing Taper Ratio
CRH = 11.1;                     % Root Chord

% Getting Swet %
SwetH = SWET(SexpH);

% Getting Skin Friction Coefficient %
CTH = sigmaH * CRH;
MACexpH = MAC(CRH,CTH)
 
RNh = ReynoldsNumber(V, MACexpH);
Cf_h = CF(RNh);

% Getting Form Factor %
Khoriztail = Kairfoil(tcH, M, sweepangleH)

% Calculating f and adding to array %
fhoriztail = F(Khoriztail, Cf_h, SwetH)
% Vertical Tail

SexpV = 161;
tcV = 0.09;
sweepangleV = 43.5;             % Sweep Angle 
sigmaV = 0.8;                   % Wing Taper Ratio
CRV = 15.5;                     % Root Chord

% Getting Swet %
SwetV = SWET(SexpV);

% Getting Skin Friction Coefficient %
CTV = sigmaV * CRV;
MACexpV = MAC(CRV,CTV)
 
RNv = ReynoldsNumber(V, MACexpV);
Cf_v = CF(RNv);

% Getting Form Factor %
Kverttail = Kairfoil(tcV, M, sweepangleV)

% Calculating f and adding to array %
fverttail = F(Kverttail, Cf_v, SwetV)
% Pylons

SwetP = 117;                    % Wetted Area
tcP = 0.06;
sweepangleP = 0;                % Sweep Angle
sigmaP = 1;                     % Taper Ratio: cT/cR
chordP = 16.2;                  % Chord

% Getting Skin Friction Coefficient %
RNp = ReynoldsNumber(V, chordP);
Cf_p = CF(RNp);

% Getting Form Factor %
Kpylon = Kairfoil(tcP, M, sweepangleP)

% Calculating f and adding to array %
fpylon = F(Kpylon, Cf_p, SwetP)
% Fuselage

Lf = 105;
Df = 11;

% Calculating Swet %
SwetF = 0.8 * pi * Df * Lf;

% Getting Skin Friction Coefficient %
RNf = ReynoldsNumber(V, Lf);
Cf_f = CF(RNf)

% Getting Form Factor %
ratioF = Lf/Df
Kfuse = KFR(ratioF);           % Via Digitized Figure 11.4

% Calculating f and adding to array %
ffuselage = F(Kfuse, Cf_f, SwetF)
% Nacelles

% Swet %
SwetN = 455;

% Getting Skin Friction Coefficient %
Ln = 16.8
RNn = ReynoldsNumber(V, 16.8);
Cf_n = CF(RNn)

% Getting Form Factor %
ratioN = 5;
Knacelle = KFR(ratioN);           % Via Digitized Figure 11.4

% Calculating f and adding to array %
fnacelle = F(Knacelle, Cf_n, SwetN)
% Total Skin Friction

ftotal = fwing + fhoriztail + fverttail + fpylon + ffuselage + fnacelle;
CDP_total = ftotal./Sref
% |Induced Drag Coeff|

% Getting CL %
W = 98000;              % Aircraft Weight
q = 0.5 * rho * (V.^2);  % Dynamic Pressure
CL = W ./ (q * Sref)     % Coeff of Lift

% Getting Aspect Ratio %
ARw = (bW^2)/Sref

% Getting CDi %
e = 1;                          % Oswald Efficiency Factor
CDi = (CL.^2) / (pi * ARw * e)   % Coeff of Induced Drag
%% |Total Drag & Lift/Drag Ratio|

ProfileDrag = CDP_total .* q .* Sref
InducedDrag = CDi .* q .* Sref
CDtotal = CDP_total + CDi
TotalDrag = CDtotal .* q .* Sref
L = W
LiftToDrag = L ./ TotalDrag
%% |Plots|

hold off;
plot(V, ProfileDrag, 'g', V, InducedDrag, 'b', V, TotalDrag, 'r')
title("Drag vs Velocity")
ylabel("Drag (lbs)")
xlabel("Velocity: (^{ft}/_{s})")
legend("Profile Drag", "Induced Drag", "Total Drag")

plot(V,LiftToDrag, 'm')
title("Lift to Drag Ratio vs V")
ylabel("Lift to Drag Ratio: ^{L}/_{D}")
xlabel("Velocity (^{ft}/_{s})")
legend("^{L}/_{D}")
%% |Functions|
% Reynolds Number Function

function RN = ReynoldsNumber(Velocity, characteristicLength)
    mu = 3.025E-7;
    rho = 0.0008754;
    V = Velocity;
    Lc = characteristicLength;
    RN = (rho * V * Lc)/mu;
end
% Swet Function

function Swet = SWET(Sexp)
    Swet = 2 * 1.02 * Sexp;
end
% MAC Function

function cbar = MAC(cR, cT)
    cbar= (2/3) * (cR + cT - ((cR*cT)/(cR+cT)));
end
% CR Exposed Function

function crexp = CREXP(cR, cT, y, b)
    crexp = cR - ((cR- cT)*(2*(y/b)));
end
% Skin Friction Coefficient

function Cf = CF(RN)
    Cf = 0.455 ./ ((log10(RN)).^2.58);
end
% Form Factor for Airfoils

function K = Kairfoil(tc, Mo, sweepAngle)
    numTerm = (2-Mo.^2) * cosd(sweepAngle);
    denTerm = sqrt(1-(Mo*cosd(sweepAngle)^2));
    Z = numTerm/denTerm;
    K = 1 + (Z * tc) + (100*tc^4);
end
% Form Factor via Fineness Ratio

function K = KFR(LbyD)
    K = 1.991*LbyD^-1.024+0.9084;
    % General model Power2:
    % Coefficients (with 95% confidence bounds):
    %    a =       1.991  (1.882, 2.101)
    %    b =      -1.024  (-1.091, -0.9582)
    %    c =      0.9084  (0.8888, 0.9279)
end
% F Function

function f = F(K, Cf, Swet)
    f = K * Cf * Swet;
end
%% 
%