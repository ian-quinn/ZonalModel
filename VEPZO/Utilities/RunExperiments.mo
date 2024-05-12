within VEPZO.Utilities;
function RunExperiments
  "Function to call the simulation tasks and mark the time elapsed"
  input String modelName="VEPZO.Samples.HybriZo"
    "Model to be simulated" annotation (dialog(enable=false));
  input Integer stopTime=86400;
  input Integer interval=600;
  output Modelica.Units.SI.Time ElapsedCPUtime;

protected
  Boolean translate;
  Real tik;
  Real tok;

algorithm
  translate := translateModel(modelName);
  // Step 2. Remove any existing results in the Variable Browser of Simulation Window
  DymolaCommands.Plot.removeResults();
  // Step 3. Remove any existing plot window open on Simulation Window
  removePlots();
  
  (ms_tik,sec_tik,min_tik,hour_tik,day_tik,mon_tik,year_tik) := 
    Modelica.Utilities.System.getTime();
  tik := (ms_tik*0.001) + (sec_tik) + (min_tik*60) + (hour_tik*3600);
  //Modelica.Utilities.Streams.print("tik = " + String(tik) + " [s]");

  simulateExtendedModel(
    modelName,
    stopTime=stopTime,
    method="dassl",
    outputInterval=interval,
    resultFile="HybriZo"); 

  (ms_tok,sec_tok,min_tok,hour_tok,day_tok,mon_tok,year_tok) :=
    Modelica.Utilities.System.getTime();
  tok := (ms_tok*0.001) + (sec_tok) + (min_tok*60) + (hour_tok*3600);
  //Modelica.Utilities.Streams.print("tok = " + String(tok) + " [s]");

  ElapsedCPUtime := abs(tok - tik);
  Modelica.Utilities.Streams.print("ElapsedCPUtime = " + String(ElapsedCPUtime)
     + " s");
end RunExperiments;