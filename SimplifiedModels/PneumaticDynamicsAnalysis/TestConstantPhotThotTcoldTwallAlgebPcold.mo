within SimplifiedModels.PneumaticDynamicsAnalysis;

model TestConstantPhotThotTcoldTwallAlgebPcold
  extends SimplifiedModels.TestFullModel(redeclare model PlantModel = SimplifiedModels.PneumaticDynamicsAnalysis.ConstantPhotThotTcoldTwallAlgebPcold, wHot(y = 0), wCold(y = if time < 5 then 0 else 0.01));
equation

annotation(
    experiment(StartTime = 0, StopTime = 4000, Tolerance = 1e-6, Interval = 8));


end TestConstantPhotThotTcoldTwallAlgebPcold;
