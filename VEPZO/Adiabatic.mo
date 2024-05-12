within VEPZO;
model Adiabatic
  Modelica.Blocks.Interfaces.RealOutput T(unit = "degC") annotation(HideResult=true);
  HeatPort port annotation(HideResult=true);
equation
  T = Modelica.Units.Conversions.to_degC(port.T);
  port.Q_flow = 0;
end Adiabatic;
