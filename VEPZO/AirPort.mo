within VEPZO;
connector AirPort "Interface for quasi one-dimensional fluid flow in a piping anology"
  package Medium = Modelica.Media.Air.SimpleAir "Medium model" annotation(choicesAllMatching = true);
  flow Medium.MassFlowRate m_flow "Mass flow rate from the connection point into the component" annotation(HideResult=true);
  flow Medium.EnthalpyFlowRate H_flow "Enthalpy flow rate into the component (if m_flow > 0, H_flow = m_flow*h)" annotation(HideResult=true);
  Medium.AbsolutePressure p "Pressure in the connection point" annotation(HideResult=true);
  Medium.SpecificEnthalpy h "Specific mixture enthalpy in the connection point" annotation(HideResult=true);
  annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end AirPort;
