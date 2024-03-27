within SimplifiedModels.NormalizationBlocks;

model OutputBlock
  
  parameter Real y_des "On-design a-priori value for y variable";
  parameter Real y_norm = y_des "Normalization value";
  final parameter Real y_calc(fixed = false);
  
  Modelica.Blocks.Interfaces.RealInput y_in( start=y_des) annotation (Placement(
      visible=true,
      transformation(
        origin={-80,0},
        extent={{-20,-20},{20,20}},
        rotation=0),
      iconTransformation(
        origin={-80,0},
        extent={{-20,-20},{20,20}},
        rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput y_out annotation (Placement(
      visible=true,
      transformation(
        origin={70,0},
        extent={{-10,-10},{10,10}},
        rotation=0),
      iconTransformation(
        origin={90,0},
        extent={{-10,-10},{10,10}},
        rotation=0)));
equation
  
   y_out = (y_in - y_calc)/y_norm;
   
initial equation
   y_calc = y_in;
   
 annotation (
    Icon(graphics={  Rectangle(fillColor = {247, 197, 159}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));

end OutputBlock;
