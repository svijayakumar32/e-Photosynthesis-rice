function Arate=EPS_Drive_GRNs(input1,CO2i,PPFDi,Tempi,Gc,GT,Einput)
global GRNC;
global GRNT;
global pcfactor;
global EnzymeAct;
global Jmax;
global BFVmax;
global FIVmax;

global cATPsyn;
global CPSi;
global cNADPHsyn;
global VfactorC;
global VfactorT;
global cpsii;
GRN_data=input1*pcfactor;

CO2in = CO2i;
Liin = PPFDi;
Tpin = Tempi;
GRNC = Gc;
GRNT = GT;
%% MOD 1 %%
%Original indices%

EnzymeAct=Einput(1:27)/30; %unit change
Jmax=EnzymeAct(27); 
BFVmax=Einput(28:45); 
FIVmax=Einput(46:66); 

% Adjust indices for enzyme activities after removing V16 % Comment out when using Jmax_adj or Vcmax_adj
% global EnzymeAct;
% EnzymeAct=Einput(1:26)/30;%unit change %EnzymeAct=Einput(1:27)/30;%unit change
% global Jmax;
% Jmax=EnzymeAct(26);% originally Jmax=EnzymeAct(27);
% global BFVmax;
% BFVmax=Einput(27:44); %BFVmax=Einput(28:45);
% global FIVmax;
% FIVmax=Einput(45:65); %FIVmax=Einput(46:66);

% %Adjust indices for enzyme activities after removing V8,V10,V16 %18.05.23
% global EnzymeAct;
% EnzymeAct=Einput(1:24)/30;%unit change %EnzymeAct=Einput(1:27)/30;%unit change
% global Jmax;
% Jmax=EnzymeAct(24);% originally Jmax=EnzymeAct(27);
% global BFVmax;
% BFVmax=Einput(25:42); %BFVmax=Einput(28:45);
% global FIVmax;
% FIVmax=Einput(43:63); %FIVmax=Einput(46:66);

if GRNC==1
    cATPsyn=GRN_data(34);%1.0447;%1.01866 WY201803
    CPSi=GRN_data(35);%1.0131;% 1.0237 WY201803
    cNADPHsyn=GRN_data(37);%1.094468408;%1.0388 WY201803
    cpsii=GRN_data(36);%1.0169;% 1.0129;%WY201803
end
if GRNC==0
    cATPsyn=1;%1.01866 WY201803
    CPSi=1;% 1.0237 WY201803
    cNADPHsyn=1;%1.0388 WY201803
    cpsii=1;% 1.0129;%WY201803   
end


VfactorC=GRN_data(1:33);

Arate=EPS_Drive(Liin,CO2in,Tpin);
end