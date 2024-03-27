within SimplifiedModels.NormalizationBlocks;

model InputBlock
  parameter Real u_des "On-design a-priori value for u_out_in variable";
  parameter Real u_norm = u_des "Normalization value";
  
  Modelica.Blocks.Interfaces.RealInput u_in annotation (Placement(
      visible=true,
      transformation(
        origin={-80,0},
        extent={{-20,-20},{20,20}},
        rotation=0),
      iconTransformation(
        origin={-80,0},
        extent={{-20,-20},{20,20}},
        rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput u_out( start=u_des) annotation (Placement(
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

  u_in = (u_out - u_des)/u_norm;
  
  annotation (
    Icon(graphics={  Rectangle(fillColor = {123, 224, 173}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}, coordinateSystem(extent = {{-100, -100}, {100, 100}})));

end InputBlock;
