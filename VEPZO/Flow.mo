within VEPZO;

model Flow "Virtual pipe connecting sub-zones"
  package Medium = Modelica.Media.Air.SimpleAir "Medium model";
  parameter Integer Direction "Mark the direction of the flow";
  constant Real g(unit = "m/s2") = 9.801 "Gravitational acceleration";
  AirPort_a port_a "Inlet port by default" annotation(Placement(transformation(extent = {{-10, 80}, {10, 100}})));
  AirPort_b port_b "Outlet port by default" annotation(Placement(transformation(extent = {{-10, -80}, {10, -100}})));
  Medium.BaseProperties medium_a "Medium at port_a";
  Medium.BaseProperties medium_b "Medium at port_b";
  Medium.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
  SI.Area A = i[1].dx * i[1].dz "Section plane area";
  SI.Density d "Mean density";
  SI.Pressure dp "Pressure drop from port_a to port_b";
  SI.Velocity velocity "Velocity of the flow";
  SI.Force F_m "Momentum force";
  SI.Force F_g "Gravitational force";
  SI.Force F_v "Viscos Force";
  Real gradU "Gradient of w velocity along x-axis";
  Real gradV "Gradient of w velocity along x-axis";
  Real gradW "Gradient of u velocity along z-axis";
  // Broadcast
  ZBeeperIn i[2] "Gather information of adjacent 2 zones" annotation(Placement(visible = true, transformation(origin = {123.63, 65.774}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  FBeeperOut o annotation(Placement(visible = true, transformation(origin = {-180, 11.693}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
initial equation
  //d = 0.5 * (medium_a.d + medium_b.d);
equation
  // Connection
  medium_a.p = port_a.p;
  medium_a.h = port_a.h;
  medium_b.p = port_b.p;
  medium_b.h = port_b.h;
  //
  d = 0.5 * (medium_a.d + medium_b.d);
  dp = medium_b.p - medium_a.p;
  // Gradients
  if Direction == 0 then
    gradW = 2 * (i[2].w - i[1].w) / (i[2].dx + i[1].dx);
    gradV = 2 * (i[2].v - i[1].v) / (i[2].dx + i[1].dx);
    gradU = 0;
  else
    if Direction == 1 then
      gradU = 2 * (i[2].u - i[1].u) / (i[2].dy + i[1].dy);
      gradW = 2 * (i[2].w - i[1].w) / (i[2].dy + i[1].dy);
      gradV = 0;
    else
      gradU = 2 * (i[2].u - i[1].u) / (i[2].dz + i[1].dz);
      gradV = 2 * (i[2].v - i[1].v) / (i[2].dz + i[1].dz);
      gradW = 0;
    end if;
  end if;
  // Driving forces
  if Direction == 0 then
    F_m = -d * A * (i[1].u ^ 2 - i[2].u ^ 2);
  else
    if Direction == 1 then
      F_m = -d * A * (i[1].v ^ 2 - i[2].v ^ 2);
    else
      F_m = -d * A * (i[1].w ^ 2 - i[2].w ^ 2);
    end if;
  end if;
  //
  if Direction == 0 then
    F_v = 0.5 * (i[1].F_vx + i[2].F_vx);
  else
    if Direction == 1 then
      F_v = 0.5 * (i[1].F_vy + i[2].F_vy);
    else
      F_v = 0.5 * (i[1].F_vz + i[2].F_vz);
    end if;
  end if;
  if Direction == 2 then
    F_g = d * A * g * 0.5 * (i[1].dz + i[2].dz);
  else
    F_g = 0;
  end if;
  der(velocity) * 0.5 * (i[1].dz + i[2].dz) * d * A = (-A * dp) + F_g + F_m + F_v;
  m_flow = d * A * velocity;
  // Design direction of mass flow rate
  m_flow = port_a.m_flow;
  // Handle reverse and zero flow
  port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, port_b.h);
  // Energy, mass balance
  port_a.m_flow + port_b.m_flow = 0;
  port_a.H_flow + port_b.H_flow = 0;
  // Broadcast
  o.gradU = gradU;
  o.gradV = gradV;
  o.gradW = gradW;
  annotation(experiment(StopTime = 100, Interval = 1), Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, origin = {50, 0}, rotation = -270, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name"), Rectangle(visible = true, fillColor = {255, 255, 255}, fillPattern = FillPattern.Backward, extent = {{-50, -100}, {50, 100}}), Line(visible = true, points = {{0, 50}, {0, -50}}), Line(visible = true, origin = {-5, -45}, points = {{5, -5}, {-5, 5}}), Line(visible = true, origin = {5, -45}, points = {{-5, -5}, {5, 5}})}));
end Flow;
