within VEPZO;
connector AirPort_b "Air connector with outlined icon"
  extends VEPZO.AirPort;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics={  Rectangle(visible = true, lineColor = {128, 128, 128}, fillColor = {255, 255, 255},
            fillPattern =                                                                                                                                                                                                        FillPattern.Solid, extent = {{-100, -100}, {100, 100}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics={  Rectangle(visible = true, lineColor = {128, 128, 128}, fillColor = {255, 255, 255},
            fillPattern =                                                                                                                                                                                                        FillPattern.Solid, extent = {{-40, -40}, {40, 40}}), Text(visible = true, origin = {0, 16.667}, textColor = {64, 64, 64}, extent = {{-160, 33.333}, {40, 73.333}}, textString = "%name")}));
end AirPort_b;
