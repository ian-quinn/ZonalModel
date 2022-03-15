within VEPZO;

model Zone "Finite volume"
  package Medium = Modelica.Media.Air.SimpleAir "Medium model";
  parameter Modelica.SIunits.DynamicViscosity mu = 1.85e-3;
  parameter Modelica.SIunits.Length dx = 1;
  parameter Modelica.SIunits.Length dy = 1;
  parameter Modelica.SIunits.Length dz = 1;
  parameter Medium.Temperature T_0 = Modelica.SIunits.Conversions.from_degC(20);
  parameter Medium.AbsolutePressure p_0 = 101325;
  parameter Boolean IsSource = false;
  //parameter Boolean Is3D;
  VEPZO.AirPort_a port_x1 annotation(Placement(transformation(extent = {{-85, -10}, {-65, 10}})));
  VEPZO.AirPort_b port_x2 annotation(Placement(transformation(extent = {{65, -10}, {85, 10}})));
  VEPZO.AirPort_a port_y1 annotation(Placement(transformation(extent = {{-35, -35}, {-15, -15}})));
  VEPZO.AirPort_b port_y2 annotation(Placement(transformation(extent = {{15, 15}, {35, 35}})));
  VEPZO.AirPort_a port_z1 annotation(Placement(transformation(extent = {{-10, -85}, {10, -65}})));
  VEPZO.AirPort_b port_z2 annotation(Placement(transformation(extent = {{-10, 65}, {10, 85}})));
  VEPZO.HeatPort port_s if IsSource annotation(Placement(visible = true, transformation(origin = {90, -90}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.SIunits.Volume V = dx * dy * dz;
  Modelica.SIunits.Velocity u "Characteristic velocity on x-axis direction";
  Modelica.SIunits.Velocity v "Characteristic velocity on x-axis direction";
  Modelica.SIunits.Velocity w "Characteristic velocity on z-axis direction";
  Medium.BaseProperties medium;
  Modelica.SIunits.Energy U;
  Modelica.SIunits.Mass m;
  Modelica.SIunits.Force F_vx;
  Modelica.SIunits.Force F_vy;
  Modelica.SIunits.Force F_vz;
  // Broadcast
  VEPZO.FBeeperIn i[6] annotation(Placement(visible = true, transformation(origin = {-86.074, 86.886}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  VEPZO.ZBeeperOut o annotation(Placement(visible = true, transformation(origin = {-174.99, 14.21}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
initial equation
  medium.p = p_0;
  medium.T = T_0;
equation
  // Connection
  medium.p = port_z1.p;
  medium.p = port_z2.p;
  medium.p = port_y1.p;
  medium.p = port_y2.p;
  medium.p = port_x1.p;
  medium.p = port_x2.p;
  medium.h = port_z1.h;
  medium.h = port_z2.h;
  medium.h = port_y1.h;
  medium.h = port_y2.h;
  medium.h = port_x1.h;
  medium.h = port_x2.h;
  if IsSource then
    medium.T = port_s.T;
  end if;
  //if Is3D then
  //  medium.p = port_y1.p;
  //  medium.p = port_y2.p;
  //  medium.h = port_y1.h;
  //  medium.h = port_y2.h;
  //end if;
  // Total quantities
  m = V * medium.d;
  U = m * medium.u;
  // Mass balance
  der(m) = port_z1.m_flow + port_z2.m_flow + port_x1.m_flow + port_x2.m_flow + port_y1.m_flow + port_y2.m_flow;
  if IsSource then
    der(U) = port_z1.H_flow + port_z2.H_flow + port_x1.H_flow + port_x2.H_flow + port_y1.H_flow + port_y2.H_flow + port_s.Q_flow;
  else
    der(U) = port_z1.H_flow + port_z2.H_flow + port_x1.H_flow + port_x2.H_flow + port_y1.H_flow + port_y2.H_flow;
  end if;
  // Characteristic velocity u
  if port_x1.m_flow > 0 then
    if port_x2.m_flow > 0 then
      u * medium.d * dy * dz = port_x1.m_flow - port_x2.m_flow;
    else
      u * medium.d * dy * dz = port_x1.m_flow;
    end if;
  else
    if port_x2.m_flow > 0 then
      u * medium.d * dy * dz = -port_x2.m_flow;
    else
      u = 0;
    end if;
  end if;
  // Characteristic velocity w
  if port_z1.m_flow > 0 then
    if port_z2.m_flow > 0 then
      w * medium.d * dx * dy = port_z1.m_flow - port_z2.m_flow;
    else
      w * medium.d * dx * dy = port_z1.m_flow;
    end if;
  else
    if port_z2.m_flow > 0 then
      w * medium.d * dx * dy = -port_z2.m_flow;
    else
      w = 0;
    end if;
  end if;
  // Characteristic velocity v
  if port_y1.m_flow > 0 then
    if port_y2.m_flow > 0 then
      v * medium.d * dx * dz = port_y1.m_flow - port_y2.m_flow;
    else
      v * medium.d * dx * dz = port_y1.m_flow;
    end if;
  else
    if port_y2.m_flow > 0 then
      v * medium.d * dx * dz = -port_y2.m_flow;
    else
      v = 0;
    end if;
  end if;
  // Viscous Force
  F_vx = mu * (i[6].gradU - i[5].gradU) * dx * dz + mu * (i[4].gradU - i[3].gradU) * dx * dy;
  F_vy = mu * (i[2].gradV - i[1].gradV) * dy * dz + mu * (i[4].gradV - i[3].gradV) * dx * dy;
  F_vz = mu * (i[2].gradW - i[1].gradW) * dy * dz + mu * (i[6].gradW - i[5].gradW) * dx * dz;
  // Broadcast
  o.dx = dx;
  o.dy = dy;
  o.dz = dz;
  o.u = u;
  o.v = v;
  o.w = w;
  o.F_vx = F_vx;
  o.F_vy = F_vy;
  o.F_vz = F_vz;
  //o.T = medium.T;
  annotation(experiment(StopTime = 3600, Interval = 1), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, origin = {-25, -25}, fillColor = {255, 255, 255}, extent = {{-75, -75}, {75, 75}}), Rectangle(visible = true, origin = {25, 25}, fillColor = {255, 255, 255}, extent = {{-75, -75}, {75, 75}}), Line(visible = true, origin = {-75, 75}, points = {{25, 25}, {-25, -25}}), Line(visible = true, origin = {75, 75}, points = {{25, 25}, {-25, -25}}), Line(visible = true, origin = {75, -75}, points = {{25, 25}, {-25, -25}}), Line(visible = true, origin = {-75, -75}, points = {{25, 25}, {-25, -25}}), Text(visible = true, origin = {0, -130}, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
end Zone;
