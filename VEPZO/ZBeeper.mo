within VEPZO;

record ZBeeper "Broadcast messages by Zone"
  SI.Length dx;
  SI.Length dy;
  SI.Length dz;
  SI.Velocity u;
  SI.Velocity v;
  SI.Velocity w;
  SI.Force F_vx;
  SI.Force F_vy;
  SI.Force F_vz;
  //SI.Temperature T;
  annotation(Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end ZBeeper;
