within VEPZO;

model Bound
  Modelica.SIunits.HeatFlowRate Q_flow "Heat flow rate from port_a -> port_b";
  HeatPort port_a;
  HeatPort port_b;
equation
  Q_flow = 0.8 * Modelica.Constants.sigma * 0.3536 * (port_a.T ^ 4 - (port_b.T - 100) ^ 4) + 0.8 * Modelica.Constants.sigma * 0.6464 * (port_a.T ^ 4 - port_b.T ^ 4);
  port_a.Q_flow = Q_flow;
  port_b.Q_flow = -Q_flow;
end Bound;
