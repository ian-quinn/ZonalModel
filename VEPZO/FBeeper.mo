within VEPZO;
record FBeeper "Broadcast messages by Flow"
  Real gradU annotation(HideResult=true);
  Real gradV annotation(HideResult=true);
  Real gradW annotation(HideResult=true);
  annotation(Diagram(coordinateSystem(extent = {{-150, -90}, {150, 90}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end FBeeper;
