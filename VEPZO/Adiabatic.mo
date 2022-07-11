within VEPZO;

model Adiabatic
  Modelica.Blocks.Interfaces.RealOutput T(unit = "degC");
  HeatPort port;
equation
  T = Modelica.SIunits.Conversions.to_degC(port.T);
  port.Q_flow = 0;
end Adiabatic;