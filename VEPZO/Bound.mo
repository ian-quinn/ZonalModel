within VEPZO;

model Bound
  SI.HeatFlowRate Q_flow "Heat flow rate from port_a -> port_b" annotation(HideResult=true);
  HeatPort port_a annotation(HideResult=true);
  HeatPort port_b annotation(HideResult=true);
equation
  Q_flow = 0.8 * Modelica.Constants.sigma * 0.3536 * (port_a.T ^ 4 - (port_b.T - 50) ^ 4) + 0.8 * Modelica.Constants.sigma * 0.6464 * (port_a.T ^ 4 - port_b.T ^ 4);
  port_a.Q_flow = Q_flow;
  port_b.Q_flow = -Q_flow;
end Bound;
