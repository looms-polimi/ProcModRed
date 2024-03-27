within SimplifiedModels.ThermalDynamicsAnalysis;

model TestConstantPhotThot
  extends SimplifiedModels.TestFullModel(redeclare model PlantModel = SimplifiedModels.ThermalDynamicsAnalysis.ConstantPhotThot, wCold(y = 0));
equation

annotation(
    experiment(StartTime = 0, StopTime = 4000, Tolerance = 1e-6, Interval = 8));
end TestConstantPhotThot;
