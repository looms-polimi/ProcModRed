within SimplifiedModels;

model TestFullModel

extends Modelica.Icons.Example;

replaceable model PlantModel = SimplifiedModels.CompleteModel;

  Modelica.Blocks.Sources.RealExpression wCold(y=if time < 5 then 0 else
0.01)
    annotation (Placement(transformation(origin = {-20, -2}, extent = {{-80, 2}, {-60, 22}})));
  Modelica.Blocks.Sources.RealExpression wHot(y=if time < 5 then 0 else 0)
    annotation (Placement(transformation(origin = {-20, -2}, extent = {{-80, -18}, {-60, 2}})));
  Modelica.Blocks.Sources.RealExpression Pel(y = completeModel.Pel) annotation (
    Placement(transformation(origin = {-10, 70}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.RealExpression TIT(y = completeModel.TIT) annotation (
    Placement(transformation(origin = {-10, 50}, extent = {{-10, -10}, {10, 10}})));
  PlantModel completeModel(pOutTurb=8000000,
    pOutHotHX=120000,
    TinCold=723.15,
    TinHot=1413.15,
    ToutColdStart=891.15,
    ToutHotStart=733.15,
    kCold=900,
    kHot=50,
    pNomCold=25000000,
    pNomHot=120000,
    wNomCold=270,
    wNomHot=62,
    gammaNomCold=3400,
    gammaNomHot=90,
    alpha=0.5,
    beta=0.8,
    cW=590,
    Mw=385000*0.6,
    Scold=6770*0.45,
    Shot=6770*0.45,
    Vcold=52,
    Vhot=13292,
    Kt=0.0045,
    eta_iso_nom=0.843,
    ToutTurbStart=763.15,
    eta_mech=0.98,
    eta_elec=0.996) annotation(
    Placement(transformation(origin = {1, -2}, extent = {{-21, -20}, {21, 20}})));
  NormalizationBlocks.InputBlock wColdNorm(u_des = 270)  annotation(
    Placement(transformation(origin = {-50, 10}, extent = {{-10, -10}, {10, 10}})));
  NormalizationBlocks.InputBlock wHotNorm(u_des = 62)  annotation(
    Placement(transformation(origin = {-50, -10}, extent = {{-10, -10}, {10, 10}})));
  NormalizationBlocks.OutputBlock powerNorm(y_des = 3.967977753682437e7)  annotation(
    Placement(transformation(origin = {30, 70}, extent = {{-10, -10}, {10, 10}})));
  NormalizationBlocks.OutputBlock TITnorm(y_des = 889.8120118760038)  annotation(
    Placement(transformation(origin = {30, 50}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(wCold.y, wColdNorm.u_in) annotation(
    Line(points = {{-78, 10}, {-58, 10}}, color = {0, 0, 127}));
  connect(wColdNorm.u_out, completeModel.wInCold) annotation(
    Line(points = {{-40, 10}, {-16, 10}}, color = {0, 0, 127}));
  connect(wHotNorm.u_out, completeModel.wInHot) annotation(
    Line(points = {{-40, -10}, {-16, -10}}, color = {0, 0, 127}));
  connect(wHot.y, wHotNorm.u_in) annotation(
    Line(points = {{-78, -10}, {-58, -10}}, color = {0, 0, 127}));
  connect(Pel.y, powerNorm.y_in) annotation(
    Line(points = {{2, 70}, {22, 70}}, color = {0, 0, 127}));
  connect(TIT.y, TITnorm.y_in) annotation(
    Line(points = {{2, 50}, {22, 50}}, color = {0, 0, 127}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    experiment(
      StopTime=4000,
      __Dymola_NumberOfIntervals=5000,
      __Dymola_Algorithm="Dassl"));

end TestFullModel;
