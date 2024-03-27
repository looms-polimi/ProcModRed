within SimplifiedModels.ThermalDynamicsAnalysis;

model ConstantPhotAlgebPcoldThotTcold
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
  Modelica.Blocks.Interfaces.RealInput wInCold annotation(
    Placement(transformation(extent = {{-120, 30}, {-80, 70}}), iconTransformation(extent = {{-100, 40}, {-60, 80}})));
  Modelica.Blocks.Interfaces.RealInput wInHot annotation(
    Placement(transformation(extent = {{-120, -70}, {-80, -30}}), iconTransformation(extent = {{-100, -60}, {-60, -20}})));
  output Types.Power Pel;
  output Types.Temperature TIT;
equation
  wInCold = wCold[1];
  wInHot = wHot[1];
  Tcold[1] = TinCold;
  Thot[1] = TinHot;
  pOutHot = pOutHotHX;
  pInCold - pCold = kCold/2*wCold[1];
  pCold - pOutCold = kCold/2*wCold[n + 1];
  pInHot - pHot = kHot/2*wNomHot;
  pHot - pOutHot = kHot/2*wNomHot;
  for i in 1:n loop
    0 = wCold[i] - wCold[i + 1];
    0 = wHot[i] - wHot[i + 1];
    0 = wCold[i + 1]*cpCold*(Tcold[i] - Tcold[i + 1]) + gammaNomCold*(pOutCold/pNomCold)^alpha*(wCold[i + 1]/wNomCold)^beta*Scold/n*(Twall[i] - (Tcold[i] + Tcold[i + 1])/2);
    0 = wHot[i + 1]*cpHot*(Thot[i] - Thot[i + 1]) + gammaNomHot*(pOutHot/pNomHot)^alpha*(wHot[i + 1]/wNomHot)^beta*Shot/n*(Twall[n + 1 - i] - (Thot[i] + Thot[i + 1])/2);
    Mw/n*cW*der(Twall[i]) = -gammaNomCold*(pOutCold/pNomCold)^alpha*(wCold[i + 1]/wNomCold)^beta*Scold/n*(Twall[i] - (Tcold[i] + Tcold[i + 1])/2) - gammaNomHot*(pOutHot/pNomHot)^alpha*(wHot[n + 1 - i]/wNomHot)^beta*Shot/n*(Twall[i] - (Thot[n + 2 - i] + Thot[n + 1 - i])/2);
  end for;
// Turbine
  wCold[n + 1] = Kt*pOutCold*sqrt(1/Rcold/Tcold[n + 1])*sqrt(1 - (1/(pOutCold/pOutTurb))^2);
  Pel = wCold[n + 1]*eta_iso_nom*cpCold*Tcold[n + 1]*(1 - (pOutTurb/pOutCold)^(Rcold/cpCold))*eta_mech*eta_elec;
  TIT = Tcold[n + 1];
initial equation
  der(Twall[:]) = zeros(n);
  annotation(
    Icon(coordinateSystem(preserveAspectRatio = false), graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {28, 108, 200}, lineThickness = 1, fillColor = {113, 113, 113}, fillPattern = FillPattern.Solid)}),
    Diagram(coordinateSystem(preserveAspectRatio = false)),
    experiment(__Dymola_NumberOfIntervals = 5000, __Dymola_Algorithm = "Dassl"));


end ConstantPhotAlgebPcoldThotTcold;
