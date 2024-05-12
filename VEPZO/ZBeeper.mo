within VEPZO;
record ZBeeper "Broadcast messages by Zone"
  SI.Length dx annotation(HideResult=true);
  SI.Length dy annotation(HideResult=true);
  SI.Length dz annotation(HideResult=true);
  SI.Velocity u annotation(HideResult=true);
  SI.Velocity v annotation(HideResult=true);
  SI.Velocity w annotation(HideResult=true);
  SI.Force F_vx annotation(HideResult=true);
  SI.Force F_vy annotation(HideResult=true);
  SI.Force F_vz annotation(HideResult=true);
  //SI.Temperature T;
  annotation(Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end ZBeeper;
