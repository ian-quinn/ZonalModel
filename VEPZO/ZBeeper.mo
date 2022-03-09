within VEPZO;

record ZBeeper "Broadcast messages by Zone"
  Modelica.SIunits.Length dx;
  Modelica.SIunits.Length dy;
  Modelica.SIunits.Length dz;
  Modelica.SIunits.Velocity u;
  Modelica.SIunits.Velocity v;
  Modelica.SIunits.Velocity w;
  Modelica.SIunits.Force F_vx;
  Modelica.SIunits.Force F_vy;
  Modelica.SIunits.Force F_vz;
  //SI.Temperature T;
  annotation(Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end ZBeeper;
