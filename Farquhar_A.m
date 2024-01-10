%% Calculate assimilation rates using curve fitting parameters in Farquhar model
Vcmax_m=105.1452508;% potato Solara Vcmax 
Jmax_m=158.6924474;% potato Solara Jmax 
Lii=1800;%light intensity
I2=Lii/2*0.85*(1-0.15);
Theta=0.7;
Gr=41.36971778;%µbar %Gamma Star from fitting
%Gr=38.6;%µbar von caemmerer 2020 %Gamma Star
Rd=2.319107419; %Rd from fitting
%Rd=1; %default value for day respiration (in light)
Kc=272;%µbar
Ko=166;%mbar Sharkey 2007
O=210;%mbar
J=(I2+Jmax_m-sqrt((I2+Jmax_m)^2-4*Theta*I2*Jmax_m))/(2*Theta);
Ci = (140:20:420)'; %Range of Ci values

% Preallocate output variables for assimilation rates
ACI_m_Vcmax=zeros(15,1); %A under Rubisco limitation
ACI_m_J=zeros(15,1); % A under electron-transport limitation
final_ACI_m=zeros(15,1); %minimum of these two A values 

% Calculate assimilation rates
for i=1:15 % No. of Ci values
ACI_m_Vcmax(i)=Vcmax_m*(Ci(i,1)-Gr)/(Ci(i,1)+Kc*(1+O/Ko))-Rd;
ACI_m_J(i)=J*(Ci(i,1)-Gr)/(4*Ci(i,1)+8*Gr)-Rd;
final_ACI_m=min(ACI_m_J,ACI_m_Vcmax);
end