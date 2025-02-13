%function LeafA=Leaf(WeatherRH,WeatherTemperature,Air_CO2,WeatherWind,Radiation_PAR,Radiation_NIR,Radiation_LW,PhotosynthesisType,Vcmax25,Jmax25,GRNC,Einput,Eio)%try including Rd as an input here
%function LeafA=Leaf(WeatherRH,WeatherTemperature,Air_CO2,WeatherWind,Radiation_PAR,Radiation_NIR,Radiation_LW,PhotosynthesisType,Vcmax25,Jmax25,GRNC,Einput,Eio,Rd,Gr)
% LeafAmatrix outputs a matrix with all convergence loop iteration values for Ci, NetAssimilation, Gs, LeafTemperature and Transpiration
function LeafAmatrix=Leaf(WeatherRH,WeatherTemperature,Air_CO2,WeatherWind,Radiation_PAR,Radiation_NIR,Radiation_LW,PhotosynthesisType,Vcmax25,Jmax25,GRNC,Einput,Eio,Rd,Gr)

PhotosynthesQ10=0;
R=8.314472E-3;%Gas constant KJ mole^{-1} K^{-1}
Convert=1E6/(2.35E5); %Convert W m^{-2} to u moles m^{-2} s^{-1}
Boltzman=5.6697E-8; % Stefan-Boltzmann constant W m^{-2} K^{-4}
LatentHeatVaporization=44000.0;%J mole^{-1}
Pressure=101325.0; % Standard atmospheric pressure Pa 
ConstantsCp=29.3;
PhotosynthesisTheta=0.76;
Rd25=1; % WY 202210 - should I modify this?
BallBerryIntercept=0.008; % can be modified to measured values
BallBerrySlope=10.5; % can be modified to measured values
%Air_CO2=400.0;
Air_O2=210.0;
WaterStressFunction=3;
WaterStressFactor=1.0;
    
MaxError = 0.25; MinError = 0.001;
ErrorCount = 0; MaxErrorCount = 1;
Previous2Gs = 0; %Stomatal conductance moles/m2 leaf area/s
Relax = 0.0; % Relaxation value for oscillation
%Initialize variables %Previous2Gs = 0; %Oscillation
[PreviousLeafState_Ci,Previous2Gs,PreviousLeafState_Gs,PreviousLeafState_Temperature,PreviousLeafMassFlux_GrossAssimilation,PreviousLeafMassFlux_NetAssimilation]=deal(zeros(1,1));
[ErrorLeafState_Ci,ErrorLeafState_Gs,ErrorLeafState_Temperature]=deal(ones(1,1));
LeafTemperature= WeatherTemperature;% Initial Leaf temperature C
Ci = 0.7 * Air_CO2;%Initial Ci u moles/mole
Gs = 0.01; % Initial stomatal conductance moles/m2 leaf area/s
Gb = 10.2; % Initial boundary layer conductance moles/m2 leaf area/s

%Make LeafA into a result matrix: dimensions depend on how many loop iterations we expect
%until convergence? 5 x ?
%LeafAmatrix=zeros(maxiterCount,5); %replace maxiterCount with max no of loop iterations
LeafAmatrix=zeros(5,5); %temporary guess for dimensions of matrix

    %Convergence loop for leaf
    max_iters = 2; % Set a limit for the maximum number of iterations to ensure we don't get stuck
    iterCount = 0; % Count number of while loop iterations
    while ((abs(ErrorLeafState_Ci) >= MinError || abs(ErrorLeafState_Gs) >= MinError ||abs(ErrorLeafState_Temperature) >= MinError) && ErrorCount <= MaxErrorCount) && iterCount < max_iters
	iterCount=iterCount+1;
            % Adding measured Rd, Gr as inputs to ComputPhotosynthesisRate
            PhotosynthesisRate=ComputPhotosynthesisRate(PhotosynthesisType,PhotosynthesQ10,Vcmax25,Jmax25,Rd25,R,LeafTemperature,Convert, Radiation_PAR,PhotosynthesisTheta,Ci,Air_O2,GRNC,Einput,Eio,Rd,Gr);
            NetAssimilation=PhotosynthesisRate(1); % move this out and calculate in adj scripts?
            GrossAssimilation=PhotosynthesisRate(2);
            Rd=PhotosynthesisRate(3);
            %GammaStar=PhotosynthesisRate(4);
            Gr=PhotosynthesisRate(4);

        %Comment out Boundary Layer Conductance calculation since it needs NetAssimilation?
        BoundaryCound=ComputeBoundaryLayerConductance(WeatherTemperature,LeafTemperature,WeatherRH,WeatherWind,Pressure,Gs,NetAssimilation,Air_CO2);
        Gb=BoundaryCound(1);
        Cb=BoundaryCound(2);
        Eb=BoundaryCound(3);
        
        % Comment out Stomatal Conductance calculation since it needs
        % NetAssimilation and inputs from Boundary Layer Conductance?
        % Also adjusts Ci so it is not exactly 0.7*Air_CO2 each time
        CalGs=ComputGsBallBerry(BallBerryIntercept,BallBerrySlope,WaterStressFactor,Eb,Cb,Gb,NetAssimilation,LeafTemperature,Air_CO2);
        Gs=CalGs(1);
        Ci=CalGs(2);
        
        % Compute energy balance
        CalLeafTemperature=ComputeEnergyBalance(Gb,Gs,NetAssimilation,LeafTemperature,WeatherRH,WeatherTemperature,Pressure,Boltzman,ConstantsCp,LatentHeatVaporization,Radiation_PAR,Radiation_NIR,Radiation_LW);
        LeafTemperature=CalLeafTemperature(1);
        Transpiration=CalLeafTemperature(2);
        
        % Compute convergence error
        ErrorLeafState_Ci = (PreviousLeafState_Ci - Ci) / Ci;
        ErrorLeafState_Gs = (PreviousLeafState_Gs - Gs) / Gs;
        ErrorLeafState_Temperature = (PreviousLeafState_Temperature - LeafTemperature) / LeafTemperature;

        %Check for oscillation and divergence to apply relaxation
        if (abs(Previous2Gs - Gs) < 0.01 && abs(ErrorLeafState_Gs) >= 0.001 && ErrorCount > 1) ||(abs(ErrorLeafState_Ci) > MaxError && ErrorCount > 1) || (abs(ErrorLeafState_Gs) > MaxError && ErrorCount > 1) ||(abs(ErrorLeafState_Temperature) > MaxError && ErrorCount > 1) % Divergence
            Relax = 0.5;
            % fprintf(LogOutputFile,"Relax");
        else if (Relax > 0.0)
            Relax = 0.0;
            end
        end
        %Apply relaxation due to oscillation and divergence
        Ci = Ci - Relax * (Ci - PreviousLeafState_Ci);
        % if (Ci <= GammaStar) % comment this loop out if we want to avoid incorporating GammaStar from ComputPhotosynthesisRate
        %     Ci = GammaStar; % ppm
        % end
        if (Ci <= Gr) % modify this loop to use Gr from A/Cc fit of measured data instead of GammaStar from ComputPhotosynthesisRate
            Ci = Gr; % ppm
        end
        Gs = Gs - Relax * (Gs - PreviousLeafState_Gs);
        if Gs < BallBerryIntercept
           Gs = BallBerryIntercept;
        end
        GrossAssimilation = GrossAssimilation - Relax *(GrossAssimilation- PreviousLeafMassFlux_GrossAssimilation);
        if (GrossAssimilation < 0.0)
            GrossAssimilation = 0.0;
            NetAssimilation = GrossAssimilation - Rd; %Comment this out if we don't calculate NetAssimilation here
        end
        %Populate matrix with results for each loop iteration (row)
        LeafAmatrix(iterCount,1)=Ci;
        LeafAmatrix(iterCount,2)=NetAssimilation;
        LeafAmatrix(iterCount,3)=Gs;
        LeafAmatrix(iterCount,4)=LeafTemperature;
        LeafAmatrix(iterCount,5)=Transpiration;
        %Update values
        PreviousLeafState_Ci = Ci;
        Previous2Gs = PreviousLeafState_Gs; %Oscillation
        PreviousLeafState_Gs = Gs;
        PreviousLeafState_Temperature = LeafTemperature;
        PreviousLeafMassFlux_GrossAssimilation = GrossAssimilation;
        PreviousLeafMassFlux_NetAssimilation = NetAssimilation; %Comment out if we don't calculate NetAssimilation here?
        ErrorCount = ErrorCount + 1;
    end

LeafA(1)=Ci;
LeafA(2)=NetAssimilation;
LeafA(3)=Gs;
LeafA(4)=LeafTemperature;
LeafA(5)=Transpiration;
end
