within SimplifiedModels.PneumaticDynamicsAnalysis;

model ConstantPhotThotTcold
parameter Integer n = 6;
  parameter Types.Pressure pOutTurb;
  parameter Types.Pressure pOutHotHX;
  parameter Types.Temperature TinCold;
  parameter Types.Temperature TinHot;
  parameter Types.Temperature ToutColdStart;
  parameter Types.Temperature ToutHotStart;
  parameter Types.LinPdropCoeff kCold;
  parameter Types.LinPdropCoeff kHot;
  parameter Types.Pressure pNomCold;
  parameter Types.Pressure pNomHot;
  parameter Types.MassFlowRate wNomCold;
  parameter Types.MassFlowRate wNomHot;
  parameter Types.CoefficientOfHeatTransfer gammaNomCold;
  parameter Types.CoefficientOfHeatTransfer gammaNomHot;
  parameter Real alpha;
  parameter Real beta;
  parameter Types.SpecificHeatCapacity cW;
  parameter Types.Mass Mw;
  parameter Types.Area Scold;
  parameter Types.Area Shot;
  parameter Types.Volume Vcold;
  parameter Types.Volume Vhot;
  parameter Types.Area Kt;
  parameter Types.Efficiency eta_iso_nom;
  parameter Types.Temperature ToutTurbStart;
  parameter Types.Efficiency eta_mech "Mechanical efficiency";
  parameter Types.Efficiency eta_elec "Nominal Electrical efficiency";
  constant Types.MolarMass MMco2 = 0.0440095;
  constant Types.MolarMass MMflue = 0.029076969929450615;
  constant Types.MassGasConstant Rcold = Modelica.Constants.R/MMco2;
  constant Types.MassGasConstant Rhot = Modelica.Constants.R/MMflue;
  constant Types.SpecificHeatCapacityAtConstantPressure cpCold = 1170;
  constant Types.SpecificHeatCapacityAtConstantPressure cpHot = 1260;
  constant Types.SpecificHeatCapacityAtConstantVolume cvCold = cpCold - Rcold;
  constant Types.SpecificHeatCapacityAtConstantVolume cvHot = cpHot - Rhot;

  Types.Pressure pCold;
  Types.Pressure pHot;
  Types.Pressure pInCold;
  Types.Pressure pInHot;
  Types.Pressure pOutCold;
  Types.Pressure pOutHot;
  Types.MassFlowRate wCold[n + 1];
  Types.MassFlowRate wHot[n + 1];
  Types.Temperature Tcold[n + 1];
  Types.Temperature Thot[n + 1];
  Types.Temperature Twall[n];

  Modelica.Blocks.Interfaces.RealInput wInCold annotation (
    Placement(transformation(extent = {{-120, 30}, {-80, 70}}), iconTransformation(extent = {{-100, 40}, {-60, 80}})));
  Modelica.Blocks.Interfaces.RealInput wInHot annotation (
    Placement(transformation(extent = {{-120, -70}, {-80, -30}}), iconTransformation(extent = {{-100, -60}, {-60, -20}})));

  output Types.Power Pel;
  output Types.Temperature TIT;
  
  Types.Time tauPrimeCold = Mw*cW/wCold[1]/cpCold;
  Types.Time tauPrimeHot = Mw*cW/wHot[1]/cpHot;

equation
  wInCold = wCold[1];
  wInHot = wHot[1];
  Tcold[1] = TinCold;
  pOutHot = pOutHotHX;
  pInCold - pCold = kCold/2*wCold[1];
  pCold - pOutCold = kCold/2*wCold[n + 1];
  pInHot - pHot = kHot/2*wNomHot;
  pHot - pOutHot = kHot/2*wNomHot;
  Thot[:] = {1140 + 273.15, 899.55786 + 273.15, 736.8236 + 273.15, 628.93317 + 273.15, 558.8643 + 273.15, 514.08795 + 273.15, 485.7986 + 273.15};
  Twall[:] = {454.55356+273.15,463.9951+273.15,478.93118+273.15,502.29544+273.15,538.2851+273.15,592.68866+273.15};

  for i in 1:n loop

    Vcold/n/(Rcold*Tcold[i+1])*der(pCold) = wCold[i] - wCold[i + 1];

    0 = wHot[i] - wHot[i + 1];
        
    pCold/Rcold/Tcold[i+1]*Vcold/n*cvCold*der(Tcold[i+1]) = wCold[i+1]*cpCold*(Tcold[i]-Tcold[i+1]) + gammaNomCold*(pOutCold/pNomCold)^alpha*(wCold[i+1]/wNomCold)^beta*Scold/n*(Twall[i] - (Tcold[i]+Tcold[i+1])/2);


  end for;
// Turbine
  wCold[n + 1] = Kt*pOutCold*sqrt(1/Rcold/Tcold[n+1])*sqrt(1 - (1/(pOutCold/pOutTurb))^2);
  Pel = wCold[n + 1] * eta_iso_nom * cpCold * Tcold[n + 1]*(1-(pOutTurb/pOutCold)^(Rcold/cpCold)) * eta_mech * eta_elec;
  TIT = Tcold[n + 1];

initial equation
  der(Tcold[:]) = zeros(n + 1);
  der(pCold) = 0;
  annotation (
    Icon(coordinateSystem(preserveAspectRatio = false), graphics={  Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {28, 108, 200}, lineThickness = 1, fillColor = {113, 113, 113}, fillPattern = FillPattern.Solid)}),
    Diagram(coordinateSystem(preserveAspectRatio = false)),
    experiment(__Dymola_NumberOfIntervals = 5000, __Dymola_Algorithm = "Dassl"));


end ConstantPhotThotTcold;
