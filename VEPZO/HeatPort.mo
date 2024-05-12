within VEPZO;
connector HeatPort "Thermal port for 1-dim. heat transfer"
  SI.Temperature T "Port temperature" annotation(HideResult=true);
  flow SI.HeatFlowRate Q_flow "Heat flow rate (positive if flowing from outside into the component)" annotation(HideResult=true);
  annotation(Documentation(info = "<html>

</html>"), Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}), graphics={  Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {191, 0, 0}, fillColor = {191, 0, 0},
            fillPattern =                                                                                                                                                                                                     FillPattern.Solid)}), Diagram(coordinateSystem(preserveAspectRatio = true, extent = {{-100, -100}, {100, 100}}, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, lineColor = {191, 0, 0}, fillColor = {191, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-50, -50}, {50, 50}}), Text(visible = true, textColor = {64, 64, 64}, extent = {{-120, 60}, {100, 120}}, textString = "%name")}));
end HeatPort;
