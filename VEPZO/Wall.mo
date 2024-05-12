within VEPZO;
model Wall
  package Medium = Modelica.Media.Air.SimpleAir "Medium model";
  parameter Integer Direction annotation(HideResult=true);
  parameter Boolean IsSource = false annotation(HideResult=true);
  parameter Boolean IsAdiabatic = true annotation(HideResult=true);
  parameter Boolean IsRadiated = false annotation(HideResult=true);
  parameter SI.Temperature T_0 = 293.15 annotation(HideResult=true);
  parameter SI.CoefficientOfHeatTransfer hc = 3.7 annotation(HideResult=true);
  parameter SI.CoefficientOfHeatTransfer hcout = 7.6 annotation(HideResult=true);
  parameter SI.ThermalResistance R = 0.1 annotation(HideResult=true);
  parameter SI.HeatCapacity C = 60000 annotation(HideResult=true);
  // cache the coordinates of the facet
  parameter SI.Length X = 0 annotation(HideResult=true);
  parameter SI.Length Y = 0 annotation(HideResult=true);
  parameter SI.Length Z = 0 annotation(HideResult=true);
  HeatPort port_s annotation(HideResult=true);
  HeatPort port_t annotation(HideResult=true);
  HeatPort port_r annotation(HideResult=true);
  //HeatPort port_s if IsSource annotation(Placement(visible = true, transformation(origin = {75, 67.532}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {80, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  //HeatPort port_t if not IsAdiabatic;
  //HeatPort port_r if IsRadiated;
  AirPort_a port_a annotation(HideResult=true, Placement(transformation(extent = {{-10, 80}, {10, 100}}, origin = {0, 0}, rotation = 0), visible = true, iconTransformation(origin = {80, 0}, extent = {{-10, 80}, {10, 100}}, rotation = 0)));
  Medium.BaseProperties medium annotation(HideResult=true);
  ZBeeperIn i annotation(HideResult=true);
  FBeeperOut o annotation(HideResult=true);
  Real gradU annotation(HideResult=true);
  Real gradV annotation(HideResult=true);
  Real gradW annotation(HideResult=true);
  SI.Area A annotation(HideResult=true);
  SI.Temperature T;
  SI.HeatFlowRate Q_flow annotation(HideResult=true);
  //SI.HeatFlowRate Q_flow if IsSource;
  //SI.Temperature Tw if IsSource;
initial equation
  T = T_0;
equation
  if Direction == 0 then
    gradV = 2 * i.v / i.dx;
    gradW = 2 * i.w / i.dx;
    gradU = 0;
    A = i.dy * i.dz;
  else
    if Direction == 1 then
      gradU = 2 * i.u / i.dy;
      gradW = 2 * i.w / i.dy;
      gradV = 0;
      A = i.dx * i.dz;
    else
      gradU = 2 * i.u / i.dz;
      gradV = 2 * i.v / i.dz;
      gradW = 0;
      A = i.dx * i.dy;
    end if;
  end if;
  port_a.p = medium.p;
  port_a.h = medium.h;
  port_a.m_flow = 0;

  port_r.T = T;
  // wall function
  medium.T - T = (R + 1 / (A * hc)) * port_a.H_flow;
  //medium.T - port_s.T = (R + 1 / (A * hc)) * Q_flow;
  port_s.T = T;
  port_a.H_flow + port_s.Q_flow = der(T) * C + Q_flow + port_r.Q_flow;

  T - port_t.T = R * Q_flow;
  port_t.Q_flow = Q_flow;

  o.gradU = gradU;
  o.gradV = gradV;
  o.gradW = gradW;
  annotation(Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics={  Text(visible = true, origin = {-20, 50}, textColor = {64, 64, 64}, extent = {{-80, 10}, {80, 50}}, textString = "%name"), Rectangle(visible = true, origin = {80, 80}, fillColor = {255, 255, 255},
            fillPattern =                                                                                                                                                                                                        FillPattern.Backward, extent = {{-20, -20}, {20, 20}})}));
end Wall;
