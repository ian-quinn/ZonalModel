model VEPZO
  model Zone "Finite volume"
    package Medium = Modelica.Media.Air.SimpleAir "Medium model";
    parameter Modelica.SIunits.DynamicViscosity mu = 1.85e-5;
    parameter Modelica.SIunits.Length dx = 1;
    parameter Modelica.SIunits.Length dy = 1;
    parameter Modelica.SIunits.Length dz = 1;
    parameter Medium.Temperature T_0 = Modelica.SIunits.Conversions.from_degC(20);
    parameter Medium.AbsolutePressure p_0 = 101325;
    //parameter Boolean Is3D;
    AirPort_b port_x1 annotation(Placement(transformation(extent = {{-85, -10}, {-65, 10}})));
    AirPort_a port_x2 annotation(Placement(transformation(extent = {{65, -10}, {85, 10}})));
    //AirPort_a port_y1 annotation(Placement(transformation(extent = {{-35, -35}, {-15, -15}})));
    //AirPort_b port_y2 annotation(Placement(transformation(extent = {{15, 15}, {35, 35}})));
    AirPort_b port_z1 annotation(Placement(transformation(extent = {{-10, -85}, {10, -65}})));
    AirPort_a port_z2 annotation(Placement(transformation(extent = {{-10, 65}, {10, 85}})));
    Modelica.SIunits.Volume V = dx * dy * dz;
    Modelica.SIunits.Velocity u "Characteristic velocity on x-axis direction";
    //Modelica.SIunits.Velocity v "Characteristic velocity on x-axis direction";
    Modelica.SIunits.Velocity w "Characteristic velocity on z-axis direction";
    Medium.BaseProperties medium;
    Modelica.SIunits.Energy U;
    Modelica.SIunits.Mass m;
    Modelica.SIunits.Force F_vx;
    //Modelica.SIunits.Force F_vy if Is3D;
    Modelica.SIunits.Force F_vz;
    // Broadcast
    FBeeperIn i[4] annotation(Placement(visible = true, transformation(origin = {-86.074, 86.886}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ZBeeperOut o annotation(Placement(visible = true, transformation(origin = {-174.99, 14.21}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  initial equation
    medium.p = p_0;
    medium.T = T_0;
  equation
    // Connection
    medium.p = port_z1.p;
    medium.p = port_z2.p;
    medium.p = port_x1.p;
    medium.p = port_x2.p;
    medium.h = port_z1.h;
    medium.h = port_z2.h;
    medium.h = port_x1.h;
    medium.h = port_x2.h;
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
    der(m) = port_z1.m_flow + port_z2.m_flow + port_x1.m_flow + port_x2.m_flow;
    der(U) = port_z1.H_flow + port_z2.H_flow + port_x1.H_flow + port_x2.H_flow;
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
    // Viscous Force
    F_vx = mu * (i[2].gradX - i[1].gradX) * dx;
    F_vz = mu * (i[4].gradZ - i[3].gradZ) * dz;
    // Broadcast
    o.dx = dx;
    o.dy = dy;
    o.dz = dz;
    o.u = u;
    o.w = w;
    o.F_vx = F_vx;
    o.F_vz = F_vz;
    annotation(experiment(StopTime = 3600, Interval = 1), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, origin = {-25, -25}, fillColor = {255, 255, 255}, extent = {{-75, -75}, {75, 75}}), Rectangle(visible = true, origin = {25, 25}, fillColor = {255, 255, 255}, extent = {{-75, -75}, {75, 75}}), Line(visible = true, origin = {-75, 75}, points = {{25, 25}, {-25, -25}}), Line(visible = true, origin = {75, 75}, points = {{25, 25}, {-25, -25}}), Line(visible = true, origin = {75, -75}, points = {{25, 25}, {-25, -25}}), Line(visible = true, origin = {-75, -75}, points = {{25, 25}, {-25, -25}}), Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
  end Zone;

  model Flow "Virtual pipe connecting sub-zones"
    package Medium = Modelica.Media.Air.SimpleAir "Medium model";
    parameter Integer Direction "Mark the direction of the flow";
    constant Real g(unit = "m/s2") = 9.801 "Gravitational acceleration";
    AirPort_a port_a "Inlet port by default" annotation(Placement(transformation(extent = {{-10, 80}, {10, 100}})));
    AirPort_b port_b "Outlet port by default" annotation(Placement(transformation(extent = {{-10, -80}, {10, -100}})));
    Medium.BaseProperties medium_a "Medium at port_a";
    Medium.BaseProperties medium_b "Medium at port_b";
    Medium.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
    Modelica.SIunits.Area A = i[1].dx * i[1].dz "Section plane area";
    Modelica.SIunits.Density d "Mean density";
    Modelica.SIunits.Pressure dp "Pressure drop from port_a to port_b";
    Modelica.SIunits.Velocity velocity "Velocity of the flow";
    Modelica.SIunits.Force F_m "Momentum force";
    Modelica.SIunits.Force F_g "Gravitational force";
    Modelica.SIunits.Force F_v "Viscos Force";
    Real gradX "Gradient of w velocity along x-axis";
    Real gradZ "Gradient of u velocity along z-axis";
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
      gradX = 2 * (i[2].w - i[1].w) / (i[2].dx + i[1].dx);
      gradZ = 0;
    else
      gradZ = 2 * (i[2].u - i[1].u) / (i[2].dz + i[1].dz);
      gradX = 0;
    end if;
    // Driving forces
    if Direction == 0 then
      F_m = -d * A * (i[1].u ^ 2 - i[2].u ^ 2);
    else
      F_m = -d * A * (i[1].w ^ 2 - i[2].w ^ 2);
    end if;
    //
    if Direction == 0 then
      F_v = 0.5 * (i[1].F_vx + i[2].F_vx);
    else
      F_v = 0.5 * (i[1].F_vz + i[2].F_vz);
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
    o.gradX = gradX;
    o.gradZ = gradZ;
    annotation(experiment(StopTime = 100, Interval = 1), Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name"), Rectangle(visible = true, fillColor = {255, 255, 255}, fillPattern = FillPattern.Backward, extent = {{-50, -100}, {50, 100}})}));
  end Flow;
public
  model Wall
    package Medium = Modelica.Media.Air.SimpleAir "Medium model";
    parameter Integer Direction;
    AirPort_a port_a annotation(Placement(transformation(extent = {{-10, 80}, {10, 100}})));
    Medium.BaseProperties medium;
    ZBeeperIn i;
    FBeeperOut o;
    Real gradX;
    Real gradZ;
  equation
    port_a.m_flow = 0;
    port_a.H_flow = 0;
    port_a.p = medium.p;
    port_a.h = medium.h;
    if Direction == 0 then
      gradX = 2 * i.w / i.dx;
      gradZ = 0;
    else
      gradZ = 2 * i.u / i.dz;
      gradX = 0;
    end if;
    o.gradX = gradX;
    o.gradZ = gradZ;
    annotation(Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 10}, {150, 50}}, textString = "%name"), Rectangle(visible = true, origin = {0, 90}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Backward, extent = {{-100, -10}, {100, 10}})}));
  end Wall;
public
  connector AirPort "Interface for quasi one-dimensional fluid flow in a piping anology"
    package Medium = Modelica.Media.Air.SimpleAir "Medium model" annotation(choicesAllMatching = true);
    flow Medium.MassFlowRate m_flow "Mass flow rate from the connection point into the component";
    flow Medium.EnthalpyFlowRate H_flow "Enthalpy flow rate into the component (if m_flow > 0, H_flow = m_flow*h)";
    Medium.AbsolutePressure p "Pressure in the connection point";
    Medium.SpecificEnthalpy h "Specific mixture enthalpy in the connection point";
    annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
  end AirPort;

  connector AirPort_a "Air connector with filled icon"
    extends AirPort;
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {128, 128, 128}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {100, 100}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {128, 128, 128}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid, extent = {{-40, -40}, {40, 40}}), Text(visible = true, origin = {0, 16.667}, textColor = {64, 64, 64}, extent = {{-160, 33.333}, {40, 73.333}}, textString = "%name")}));
  end AirPort_a;

  connector AirPort_b "Air connector with outlined icon"
    extends AirPort;
    annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {128, 128, 128}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-100, -100}, {100, 100}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {128, 128, 128}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-40, -40}, {40, 40}}), Text(visible = true, origin = {0, 16.667}, textColor = {64, 64, 64}, extent = {{-160, 33.333}, {40, 73.333}}, textString = "%name")}));
  end AirPort_b;

  record ZBeeper "Broadcast messages by Zone"
    Modelica.SIunits.Length dx;
    Modelica.SIunits.Length dy;
    Modelica.SIunits.Length dz;
    Modelica.SIunits.Velocity u;
    Modelica.SIunits.Velocity w;
    Modelica.SIunits.Force F_vx;
    Modelica.SIunits.Force F_vz;
    annotation(Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
  end ZBeeper;

  connector ZBeeperIn = input ZBeeper annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Ellipse(visible = true, fillColor = {255, 255, 255}, extent = {{-50, -50}, {50, 50}}), Text(visible = true, extent = {{-40, -40}, {40, 40}}, textString = "Z")}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {128, 128, 128}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid, extent = {{-40, -40}, {40, 40}}), Text(visible = true, origin = {0, 16.667}, textColor = {64, 64, 64}, extent = {{-160, 33.333}, {40, 73.333}}, textString = "%name")}));
  connector ZBeeperOut = output ZBeeper annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Ellipse(visible = true, fillPattern = FillPattern.Solid, extent = {{-50, -50}, {50, 50}}), Text(visible = true, textColor = {255, 255, 255}, extent = {{-40, -40}, {40, 40}}, textString = "Z")}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {128, 128, 128}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid, extent = {{-40, -40}, {40, 40}}), Text(visible = true, origin = {0, 16.667}, textColor = {64, 64, 64}, extent = {{-160, 33.333}, {40, 73.333}}, textString = "%name")}));

  record FBeeper "Broadcast messages by Flow"
    Real gradX;
    Real gradZ;
    annotation(Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
  end FBeeper;

  connector FBeeperIn = input FBeeper annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Ellipse(visible = true, fillColor = {255, 255, 255}, extent = {{-50, -50}, {50, 50}}), Text(visible = true, extent = {{-40, -40}, {40, 40}}, textString = "F")}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {128, 128, 128}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid, extent = {{-40, -40}, {40, 40}})}));
  connector FBeeperOut = output FBeeper annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Ellipse(visible = true, fillPattern = FillPattern.Solid, extent = {{-50, -50}, {50, 50}}), Text(visible = true, textColor = {255, 255, 255}, extent = {{-40, -40}, {40, 40}}, textString = "F")}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {128, 128, 128}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid, extent = {{-40, -40}, {40, 40}})}));
public
  model Test
    Zone zone1(T_0 = Modelica.SIunits.Conversions.from_degC(80)) annotation(Placement(visible = true, transformation(origin = {50, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Zone zone2 annotation(Placement(visible = true, transformation(origin = {0, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Zone zone3 annotation(Placement(visible = true, transformation(origin = {-50, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Zone zone4 annotation(Placement(visible = true, transformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Zone zone5 annotation(Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Zone zone6 annotation(Placement(visible = true, transformation(origin = {-50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Zone zone7 annotation(Placement(visible = true, transformation(origin = {50, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Zone zone8 annotation(Placement(visible = true, transformation(origin = {-0, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Zone zone9 annotation(Placement(visible = true, transformation(origin = {-50, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Flow flow1(Direction = 0) annotation(Placement(visible = true, transformation(origin = {25, -50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Flow flow2(Direction = 0) annotation(Placement(visible = true, transformation(origin = {-25, -50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Flow flow3(Direction = 2) annotation(Placement(visible = true, transformation(origin = {50, -25}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Flow flow4(Direction = 2) annotation(Placement(visible = true, transformation(origin = {0, -25}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Flow flow5(Direction = 2) annotation(Placement(visible = true, transformation(origin = {-50, -25}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Flow flow6(Direction = 0) annotation(Placement(visible = true, transformation(origin = {25, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Flow flow7(Direction = 0) annotation(Placement(visible = true, transformation(origin = {-25, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Flow flow8(Direction = 2) annotation(Placement(visible = true, transformation(origin = {50, 25}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Flow flow9(Direction = 2) annotation(Placement(visible = true, transformation(origin = {0, 25}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Flow flow10(Direction = 2) annotation(Placement(visible = true, transformation(origin = {-50, 25}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Flow flow11(Direction = 0) annotation(Placement(visible = true, transformation(origin = {25, 50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Flow flow12(Direction = 0) annotation(Placement(visible = true, transformation(origin = {-25, 50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Wall wall1(Direction = 2) annotation(Placement(visible = true, transformation(origin = {50, -75}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Wall wall2(Direction = 2) annotation(Placement(visible = true, transformation(origin = {0, -75}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Wall wall3(Direction = 2) annotation(Placement(visible = true, transformation(origin = {-50, -75}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Wall wall4(Direction = 0) annotation(Placement(visible = true, transformation(origin = {-75, -50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Wall wall5(Direction = 0) annotation(Placement(visible = true, transformation(origin = {-75, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Wall wall6(Direction = 0) annotation(Placement(visible = true, transformation(origin = {-75, 50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Wall wall7(Direction = 2) annotation(Placement(visible = true, transformation(origin = {-50, 75}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
    Wall wall8(Direction = 2) annotation(Placement(visible = true, transformation(origin = {0, 75}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
    Wall wall9(Direction = 2) annotation(Placement(visible = true, transformation(origin = {50, 75}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));
    Wall wall10(Direction = 0) annotation(Placement(visible = true, transformation(origin = {75, 50}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
    Wall wall11(Direction = 0) annotation(Placement(visible = true, transformation(origin = {75, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
    Wall wall12(Direction = 0) annotation(Placement(visible = true, transformation(origin = {75, -50}, extent = {{-10, -10}, {10, 10}}, rotation = -270)));
  equation
    connect(zone1.port_x1, flow1.port_a) annotation(Line(visible = true, origin = {38.25, -50}, points = {{4.25, 0}, {-4.25, 0}}, color = {128, 128, 128}));
    connect(zone1.port_z2, flow3.port_b) annotation(Line(visible = true, origin = {50, -38.25}, points = {{0, -4.25}, {0, 4.25}}, color = {128, 128, 128}));
    connect(flow3.port_a, zone4.port_z1) annotation(Line(visible = true, origin = {50, -11.75}, points = {{-0, -4.25}, {0, 4.25}}, color = {128, 128, 128}));
    connect(zone4.port_x1, flow6.port_a) annotation(Line(visible = true, origin = {38.25, 0}, points = {{4.25, -0}, {-4.25, 0}}, color = {128, 128, 128}));
    connect(zone4.port_z2, flow8.port_b) annotation(Line(visible = true, origin = {50, 11.75}, points = {{0, -4.25}, {-0, 4.25}}, color = {128, 128, 128}));
    connect(zone7.port_z1, flow8.port_a) annotation(Line(visible = true, origin = {50, 38.25}, points = {{0, 4.25}, {0, -4.25}}, color = {128, 128, 128}));
    connect(zone7.port_x1, flow11.port_a) annotation(Line(visible = true, origin = {38.25, 50}, points = {{4.25, 0}, {-4.25, 0}}, color = {128, 128, 128}));
    connect(zone8.port_x2, flow11.port_b) annotation(Line(visible = true, origin = {11.75, 50}, points = {{-4.25, 0}, {4.25, -0}}, color = {128, 128, 128}));
    connect(zone8.port_x1, flow12.port_a) annotation(Line(visible = true, origin = {-11.75, 50}, points = {{4.25, 0}, {-4.25, 0}}, color = {128, 128, 128}));
    connect(zone8.port_z1, flow9.port_a) annotation(Line(visible = true, origin = {0, 38.25}, points = {{-0, 4.25}, {0, -4.25}}, color = {128, 128, 128}));
    connect(zone5.port_x2, flow6.port_b) annotation(Line(visible = true, origin = {11.75, 0}, points = {{-4.25, 0}, {4.25, 0}}, color = {128, 128, 128}));
    connect(zone5.port_z2, flow9.port_b) annotation(Line(visible = true, origin = {0, 11.75}, points = {{0, -4.25}, {0, 4.25}}, color = {128, 128, 128}));
    connect(zone5.port_x1, flow7.port_a) annotation(Line(visible = true, origin = {-11.75, 0}, points = {{4.25, 0}, {-4.25, -0}}, color = {128, 128, 128}));
    connect(zone5.port_z1, flow4.port_a) annotation(Line(visible = true, origin = {0, -11.75}, points = {{0, 4.25}, {0, -4.25}}, color = {128, 128, 128}));
    connect(zone2.port_x2, flow1.port_b) annotation(Line(visible = true, origin = {11.75, -50}, points = {{-4.25, 0}, {4.25, 0}}, color = {128, 128, 128}));
    connect(zone2.port_z2, flow4.port_b) annotation(Line(visible = true, origin = {0, -38.25}, points = {{0, -4.25}, {0, 4.25}}, color = {128, 128, 128}));
    connect(zone2.port_x1, flow2.port_a) annotation(Line(visible = true, origin = {-11.75, -50}, points = {{4.25, 0}, {-4.25, 0}}, color = {128, 128, 128}));
    connect(zone3.port_x2, flow2.port_b) annotation(Line(visible = true, origin = {-35.677, -50.741}, points = {{-6.823, 0.741}, {3.411, 0.741}, {1.677, 0.741}}, color = {128, 128, 128}));
    connect(zone3.port_z2, flow5.port_b) annotation(Line(visible = true, origin = {-50, -38.25}, points = {{0, -4.25}, {0, 4.25}}, color = {128, 128, 128}));
    connect(zone6.port_x2, flow7.port_b) annotation(Line(visible = true, origin = {-38.25, 0}, points = {{-4.25, 0}, {4.25, -0}}, color = {128, 128, 128}));
    connect(zone6.port_z2, flow10.port_b) annotation(Line(visible = true, origin = {-50, 11.75}, points = {{0, -4.25}, {0, 4.25}}, color = {128, 128, 128}));
    connect(zone6.port_z1, flow5.port_a) annotation(Line(visible = true, origin = {-50, -11.75}, points = {{0, 4.25}, {0, -4.25}}, color = {128, 128, 128}));
    connect(zone9.port_x2, flow12.port_b) annotation(Line(visible = true, origin = {-38.25, 50}, points = {{-4.25, 0}, {4.25, 0}}, color = {128, 128, 128}));
    connect(zone9.port_z1, flow10.port_a) annotation(Line(visible = true, origin = {-50, 38.25}, points = {{0, 4.25}, {0, -4.25}}, color = {128, 128, 128}));
    connect(zone1.port_z1, wall1.port_a) annotation(Line(visible = true, origin = {50, -61.75}, points = {{0, 4.25}, {0, -4.25}}, color = {128, 128, 128}));
    connect(zone2.port_z1, wall2.port_a) annotation(Line(visible = true, origin = {0, -61.75}, points = {{0, 4.25}, {0, -4.25}}, color = {128, 128, 128}));
    connect(zone3.port_z1, wall3.port_a) annotation(Line(visible = true, origin = {-50, -61.75}, points = {{0, 4.25}, {-0, -4.25}}, color = {128, 128, 128}));
    connect(zone3.port_x1, wall4.port_a) annotation(Line(visible = true, origin = {-63.75, -50}, points = {{6.25, 0}, {-2.25, 0}}, color = {128, 128, 128}));
    connect(zone6.port_x1, wall5.port_a) annotation(Line(visible = true, origin = {-61.75, 0}, points = {{4.25, 0}, {-4.25, 0}}, color = {128, 128, 128}));
    connect(zone9.port_x1, wall6.port_a) annotation(Line(visible = true, origin = {-61.75, 50}, points = {{4.25, 0}, {-4.25, 0}}, color = {128, 128, 128}));
    connect(zone9.port_z2, wall7.port_a) annotation(Line(visible = true, origin = {-50, 61.75}, points = {{0, -4.25}, {0, 4.25}}, color = {128, 128, 128}));
    connect(zone8.port_z2, wall8.port_a) annotation(Line(visible = true, origin = {0, 61.75}, points = {{0, -4.25}, {0, 4.25}}, color = {128, 128, 128}));
    connect(zone7.port_z2, wall9.port_a) annotation(Line(visible = true, origin = {50, 61.75}, points = {{0, -4.25}, {0, 4.25}}, color = {128, 128, 128}));
    connect(zone7.port_x2, wall10.port_a) annotation(Line(visible = true, origin = {61.75, 50}, points = {{-4.25, 0}, {4.25, 0}}, color = {128, 128, 128}));
    connect(zone4.port_x2, wall11.port_a) annotation(Line(visible = true, origin = {61.75, 0}, points = {{-4.25, 0}, {4.25, 0}}, color = {128, 128, 128}));
    connect(zone1.port_x2, wall12.port_a) annotation(Line(visible = true, origin = {62.847, -50}, points = {{-5.347, 0}, {3.153, 0}}, color = {128, 128, 128}));
    // beepers
    connect(flow1.i[1], zone2.o);
    connect(flow1.i[2], zone1.o);
    connect(flow2.i[1], zone3.o);
    connect(flow2.i[2], zone2.o);
    connect(flow3.i[1], zone1.o);
    connect(flow3.i[2], zone4.o);
    connect(flow4.i[1], zone2.o);
    connect(flow4.i[2], zone5.o);
    connect(flow5.i[1], zone3.o);
    connect(flow5.i[2], zone6.o);
    connect(flow6.i[1], zone5.o);
    connect(flow6.i[2], zone4.o);
    connect(flow7.i[1], zone6.o);
    connect(flow7.i[2], zone5.o);
    connect(flow8.i[1], zone4.o);
    connect(flow8.i[2], zone7.o);
    connect(flow9.i[1], zone5.o);
    connect(flow9.i[2], zone8.o);
    connect(flow10.i[1], zone6.o);
    connect(flow10.i[2], zone9.o);
    connect(flow11.i[1], zone8.o);
    connect(flow11.i[2], zone7.o);
    connect(flow12.i[1], zone9.o);
    connect(flow12.i[2], zone8.o);
    //
    connect(zone1.i[1], flow1.o);
    connect(zone1.i[2], wall2.o);
    connect(zone1.i[3], wall1.o);
    connect(zone1.i[4], flow3.o);
    connect(zone2.i[1], flow2.o);
    connect(zone2.i[2], flow1.o);
    connect(zone2.i[3], wall2.o);
    connect(zone2.i[4], flow4.o);
    connect(zone3.i[1], wall4.o);
    connect(zone3.i[2], flow2.o);
    connect(zone3.i[3], wall3.o);
    connect(zone3.i[4], flow5.o);
    connect(zone4.i[1], flow6.o);
    connect(zone4.i[2], wall11.o);
    connect(zone4.i[3], flow3.o);
    connect(zone4.i[4], flow8.o);
    connect(zone5.i[1], flow7.o);
    connect(zone5.i[2], flow6.o);
    connect(zone5.i[3], flow4.o);
    connect(zone5.i[4], flow9.o);
    connect(zone6.i[1], wall5.o);
    connect(zone6.i[2], flow7.o);
    connect(zone6.i[3], flow5.o);
    connect(zone6.i[4], flow10.o);
    connect(zone7.i[1], flow11.o);
    connect(zone7.i[2], wall10.o);
    connect(zone7.i[3], flow8.o);
    connect(zone7.i[4], wall9.o);
    connect(zone8.i[1], flow12.o);
    connect(zone8.i[2], flow11.o);
    connect(zone8.i[3], flow9.o);
    connect(zone8.i[4], wall8.o);
    connect(zone9.i[1], wall6.o);
    connect(zone9.i[2], flow12.o);
    connect(zone9.i[3], flow10.o);
    connect(zone9.i[4], wall7.o);
    //
    connect(wall1.i, zone1.o);
    connect(wall2.i, zone2.o);
    connect(wall3.i, zone3.o);
    connect(wall4.i, zone3.o);
    connect(wall5.i, zone6.o);
    connect(wall6.i, zone9.o);
    connect(wall7.i, zone9.o);
    connect(wall8.i, zone8.o);
    connect(wall9.i, zone7.o);
    connect(wall10.i, zone7.o);
    connect(wall11.i, zone4.o);
    connect(wall12.i, zone1.o);
    annotation(experiment(StopTime = 100, Interval = 1, __Wolfram_Algorithm = "dassl"), Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {0, 114, 195}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}, radius = 25), Text(visible = true, textColor = {64, 64, 64}, extent = {{-150, 110}, {150, 150}}, textString = "%name")}));
  end Test;
  annotation(experiment(StopTime = 3000, Interval = 1), Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {128, 0, 0}, fillColor = {255, 0, 0}, pattern = LinePattern.None, fillPattern = FillPattern.Sphere, extent = {{-100, -100}, {100, 100}}, radius = 25), Polygon(visible = true, origin = {-0, -0.523}, lineColor = {255, 255, 255}, fillColor = {255, 255, 255}, lineThickness = 7, points = {{0, 67.387}, {60, 40.523}, {60, -39.477}, {0, -69.477}, {-60, -39.477}, {-60, 40.523}}), Line(visible = true, origin = {0, -10}, points = {{-60, 10}, {0, -20}, {60, 10}}, color = {255, 255, 255}, thickness = 7), Line(visible = true, origin = {-40, 3.333}, points = {{10, -58.333}, {10, 21.667}, {-20, 36.667}}, color = {255, 255, 255}, thickness = 7), Line(visible = true, origin = {40, 3.333}, points = {{-10, -58.333}, {-10, 21.667}, {20, 36.667}}, color = {255, 255, 255}, thickness = 7), Line(visible = true, origin = {0, 40}, points = {{-30, -15}, {30, 15}}, color = {255, 255, 255}, thickness = 7), Line(visible = true, origin = {0, 40}, points = {{-30, 15}, {30, -15}}, color = {255, 255, 255}, thickness = 7), Line(visible = true, origin = {0, -50}, points = {{0, 20}, {0, -20}}, color = {255, 255, 255}, thickness = 7), Line(visible = true, origin = {0, 20}, points = {{0, 20}, {0, -20}}, color = {255, 255, 255}, thickness = 7), Line(visible = true, origin = {0, -10}, points = {{-30, -5}, {0, 10}, {30, -5}}, color = {255, 255, 255}, thickness = 7)}));
end VEPZO;
