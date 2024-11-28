%   Copyright   Xin-Guang Zhu and Stephen P. Long, University of Illinois 
%   Copyright ©  2007

%   This file is part of CarbonMetabolism.

%    CarbonMetabolism is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.

%    CarbonMetabolism is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License (GPL)
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Vrubusco_adj;
Vrubusco_adj=1.0; %keep as 1.0 after fixed optimization
global VmaxAdj;%
VmaxAdj=1.0;
Ci=[129,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420]; %CO2 input values
global pcfactor; 
global VmaxNum;

global GP;% This is a parameter used for transfer information to the PSRate file. 
% The following variables are used to output data from the simulation. 
global MW;
global BK;
global d_plot;
global Tt_plot;
global TimeBegin;
global Optiona;

global NTotal;
for CO2i = Ci %loop through each CO2 input value
    %%%%%%%%%%%%%%%%%%% Initialize Variables %%%%%%%%%%%%%%%%%%%
    % The option used in the mutate.m file for creating mutation.

    pcfactor=1;
    GP = 1;         
    TimeBegin = clock;
    Optiona =0;     

    %WARNING population size (population_Size) must be divisible by 4
    population_Size = 16; %default=16    
    number_of_generations = 1500; % The total number of generations default = 1500 %test i=100  
    mutatePercentage = 0.02; % Maximal percentage of changes in Vmax in each generation - lower to 0.0002?

    % Different options of generating new population. Default is 1. 
    generationTransfer = 1; 
    VmaxNum = 23;           % The number of enzymes used in the optimization. %originally 26 
    factor = 1;
    ScaleR = factor;        % This is a factor used to modify the enzyme concentration.

    %%%% Initialize Population Array %%%%

    % Coeff = 1/30/3600*1000;
    population= zeros(VmaxNum+2,population_Size); % shape 25xPopulation
    PPFDi=2000; %PPFD from IRRI
    WeatherTemp=28.9310407291759;
    GRNC=0;
    Einput=ones(37,1);%No gene expression data input
    %Edata=importdata('Einput7.txt');
    Edata=importdata('Einput_rice.txt');
    Eio=Edata.data(:,1);
    Eio(1)=Edata.data(1,1)*Vrubusco_adj;
    Eio(2:26)=Edata.data(2:26,1)*VmaxAdj;
    population(3:25,:) = repmat(Eio([1:6,8,10,11,13:26]),[1,population_Size]); %taking out V8, V10 and V16 in optimization
    
    % disp('Initial population:');population(3:25,:)
    %%%% Initialize BK and MW Array %%%%
    Enzyme=importdata('MW&Kcat.txt');
    MWKcat=Enzyme.data;
    MWKcat([7,9,12],:) = []; %remove rows corresponding to V8, V10 and V16
    BK=MWKcat(:,2);
    %BK([7,9,12],:) = [];
    MW=MWKcat(:,3);


    totals=population(3:VmaxNum+2,1);
    totals(1,1)=1.2*totals(1,1);
    Bks=BK(1:VmaxNum);
    Mws=MW(1:VmaxNum);
    sumd = sum(totals./Bks.*Mws);
    %%%% Mutate Initial Population %%%%
    tempn = population_Size/8;
    temp_pop=population(:,1:tempn);
    population= mutate(population, population_Size, mutatePercentage);
    population(:,1:tempn) = temp_pop;

    CO2PLOT = zeros(3,5);  
    % 1: Generation number; 2: Average; 3, Error bar; 

    TAVMatrix = zeros(number_of_generations,VmaxNum+2);
    TSTDMatrix = zeros(number_of_generations,VmaxNum+2);
    TBestMatrix = zeros(number_of_generations,VmaxNum+2);

    for i = 1:number_of_generations
        disp('generation:');i
        PopulationSamples=repmat(Eio([1:6,8,10,11,13:66]),[1,population_Size]);
        PopulationSamples(1:23,:)=population(3:25,:);
        population_idx=1:population_Size;
        mask_Glucose1PhosphateAdenylyltransferase = PopulationSamples(9,:)<Eio(11);
        mask_FructoseBisphosphateAldolaseC = PopulationSamples(17,:)<Eio(20);
        mask_FructoseBisphosphataseC = PopulationSamples(18,:)<Eio(21);
        mask_UTPglucose1PhosphateUridylyltransferase = PopulationSamples(19,:)<Eio(22);
        mask_SucrosePhosphateSynthase = PopulationSamples(20,:)<Eio(23);
        mask_SucrosePhosphatePhosphatase = PopulationSamples(21,:)<Eio(24);
        mask_Fructose26Bisphosphate2Phosphatase = PopulationSamples(22,:)<Eio(25);
        mask_6Phosphofructo2Kinase = PopulationSamples(23,:)<Eio(26);
        mask = mask_Glucose1PhosphateAdenylyltransferase | mask_FructoseBisphosphateAldolaseC | mask_FructoseBisphosphataseC | mask_UTPglucose1PhosphateUridylyltransferase | mask_SucrosePhosphateSynthase | mask_SucrosePhosphatePhosphatase | mask_Fructose26Bisphosphate2Phosphatase | mask_6Phosphofructo2Kinase;
        PopulationSamples(:,mask) = []; %remove samples that do not meet constraints
        population_idx(mask) = []; 
        size_PopulationSamples_in_dim1 = size(PopulationSamples,2);
        samples_phot = vertcat(PopulationSamples(1:6,:),PopulationSamples(4,:),PopulationSamples(7,:),PopulationSamples(6,:),PopulationSamples(8:9,:),repmat(Eio(12),[1,size_PopulationSamples_in_dim1]),PopulationSamples(10:63,:)); %remove V8, V10 and V16
        global tglobal;

        for j = 1:size_PopulationSamples_in_dim1
           
            Temp=EPS_Drive_GRNs(Einput,CO2i,PPFDi,WeatherTemp,GRNC,0,samples_phot(:,j));
            sizeT=size(d_plot);

            Tcheck = tglobal * 4/5;  
            Tindex = find(Tt_plot>Tcheck);    
            IndexMiddle = Tindex(1);
            ddiff=d_plot(sizeT(1),53:87)-d_plot(IndexMiddle,53:87);%metabolites
            dnorm1 = norm(ddiff);
            tdiff = Tt_plot(sizeT(1))-Tt_plot(IndexMiddle);
            sloped = dnorm1/tdiff/0.01;
            
            if sloped < 10^(-3)
                population(2,population_idx(j)) = Temp;
            else
                population(2,population_idx(j)) = 0;
            end
            
            
            test = population(2,population_idx(j));
   
        end
        %%%%%%%%%%% Rank Population %%%%%%%%%%%
        population= rankPop(population,population_Size);
        
        %%%%%% Save optimization Statistics %%%%%%%%
        AVMatrix = average(population,population_Size);
        STDMatrix = stdev(population,population_Size);
        BestMatrix = population(:,1)';
        
        TAVMatrix(i,:)=AVMatrix;
        TSTDMatrix(i,:)=STDMatrix;
        TBestMatrix(i,:)=BestMatrix;
        
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Check if this is the last generation
        CO2PLOT(1,i) = i;
        CO2PLOT(2,i) = AVMatrix(2);
        CO2PLOT(3,i) = STDMatrix(2);


        switch generationTransfer
            case 1
                population= resizePop(population, population_Size);
            case 2
                population= twoPoint(population,population_Size);
        end
        % Mutate Population
        population= mutate(population, population_Size, mutatePercentage);
	    if rem(i/100,1)==0
            i
            Ci_str = num2str(CO2i);
            task_id=getenv('SLURM_ARRAY_TASK_ID');
            % Specify unique filenames for each Ci
            workspacefileName = strcat ("CO2_rice_",Ci_str,"_",task_id,".mat");
            % Save the work space 
            save(workspacefileName);
            %Save matrix of optimal enzyme rates to output file
            BestMatrix=BestMatrix'; % Transpose matrix
            BestMatrixfileName = strcat ("outputenz_",Ci_str,"_",task_id,".txt");
            d_plotfileName = strcat ("d_plot_",Ci_str,"_",task_id,".xls");
            writematrix(BestMatrix,BestMatrixfileName);
            writematrix(d_plot,d_plotfileName);
        end
    end
end
%sumd = sumd- pop(9,1)/BK(7)*MW(7)-pop(11,1)/BK(9)*MW(9)-pop(13,1)/BK(11)*MW(11)-pop(14,1)/BK(12)*MW(12);      % The transketolase is double counted. So, it is corrected at here. 


% NConcDefault = sumd/1000/33 * 0.16;        % This is the default nitrogen concentration for those enzymes in the carbon metabolism pathway. 
% NConc = NConcDefault * factor;
% NConc = NConc * 1000 * 33 /0.16       % mg protein l-1

%NTotal = sumd * factor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global gp2condition_RuCon;

    
    %%%%%%%%%%%%% Test Population %%%%%%%%%%%%%%
    
% Store CO2 Assimilation Rate into pop
    
      
%     display('Time passed:'); 
%     Lapse = clock-TimeBegin
    
    
%         gp2condition_RuCon = pop(3,j)/2;
%         global gp2V111;
%         gp2V111 = pop(3,j) * 0.24;            % Here the rate of Rubisco oxygenation is adjusted.
            %Temp = CM_Drive(pop, j);
             %adjust to 23 enzymes
%           Eiopop(7)=Eiopop(4); %These duplicates are already missing from the optimization list so no need to set them as equal
%           Eiopop(9)=Eiopop(6);
% 	        Eiopop(11)=Eio(11);
%           Eiopop(12)=Eio(12);%VATP synthesis not used
            %Fixing minimal enzyme level constraints for Glucose-1-phosphate adenylyltransferase 
            % and enzymes involved in sucrose metabolism
            % from starting enzyme levels in Einput_rice
             % All values within constraints, so analysis can proceed
            %Use Eiopop from optimization loop, setting V7=V4 and V10=V7 
%             Eiopop_phot=vertcat(Eiopop(1:6),Eiopop(4),Eiopop(7),Eiopop(6),Eiopop(8:63));
%             Temp=EPS_Drive_GRNs(Einput,CO2i,PPFDi,WeatherTemp,GRNC,0,Eiopop_phot);
            % Set ATP synthase same as original
            

   

    %mcc nojvm does not use the Java Virtual Machine (JVM) so comment out
    %all figure rendering commands to avoid errors
    %errorbar(CO2PLOT(1,:),CO2PLOT(2,:),CO2PLOT(3,:)); pause(5);
    %xlim([1, i+1]); ylim([0,80]);  
        %Resize Pop Array
       