within VEPZO.Utilities;
function MatrixConverter
  "Function to import trajectory result files and write them as Matlab compatible .mat files"
  input String filename="mrt_singapore.mat" "File to be converted";
  input String outputFilename="mrt_singapore_temp.mat";
  input Integer timeSpan=600 "Time span for snapshot retrieving (seconds)"
    annotation (Dialog(__Dymola_loadSelector(filter="Matlab files (*.mat)",
    caption="Select the results trajectory file")));

protected
  Integer noRows "Number of rows in the trajectory being converted";
  Integer matrixCol "Number of columns in the filtered matrix";
  Integer matrixRow "Number of rows in the filtered matrix";
  Integer noSteps "Number of all snapshots";
  String varNames[:] "All variable names in the trajectory";
  String varCacheNames[:] "Temperal cache for interested variable names";
  String varFilteredNames[:] "Interested variable names";
  Real data[:,:] "Data read in from trajectory file";
  Real dataDump[:,:] "Sacrificial dump variable for writeMatrix command";
  Integer timeStep=1 "Time counter for snapshot retrieving";
  Integer i=1 "Loop counter";
  Integer varCount=2 "Interested variable counter";
  Integer zoneCount=0 "Record the number of zones";

algorithm
  noRows := DymolaCommands.Trajectories.readTrajectorySize(filename);
  varNames := DymolaCommands.Trajectories.readTrajectoryNames(filename);
  noColumns := size(varNames, 1);
  Modelica.Utilities.Streams.print("Length: " + String(noColumns));

  // expand the memory for cache
  varCacheNames := fill("", 5000);
  varCacheNames[1] := "Time";
  varCacheNames[2] := "ambient.y";
  i := 1;
  while i <= noColumns loop
    if Modelica.Utilities.Strings.find(varNames[i], "zone") <> 0 then
      if Modelica.Utilities.Strings.find(varNames[i], ".medium.T") <> 0 then
        if Modelica.Utilities.Strings.find(varNames[i], ".medium.T_degC") == 0 then
          // start from the 3rd column, previous two are "Time" "ambient.y"
          varCount := varCount + 1;
          varCacheNames[varCount] := varNames[i];
        end if;
      end if;
    end if;
    i := i + 1;
  end while;

  zoneCount := varCount;

  i := 1;
  while i <= noColumns loop
    if Modelica.Utilities.Strings.find(varNames[i], "wall") <> 0 then
      if Modelica.Utilities.Strings.find(varNames[i], "wallR") == 0 then
        if Modelica.Utilities.Strings.find(varNames[i], ".T") <> 0 then
          varCount := varCount + 1;
          varCacheNames[varCount] := varNames[i];
        end if;
      end if;
    end if;
    i := i + 1;
  end while;

  // note that the varCount will be "Time" plus all variable names
  // .mat file must include the retrieving variable names by readTrajectory
  varFilteredNames := varCacheNames[1:varCount];

  data := DymolaCommands.Trajectories.readTrajectory(
    filename,
    varFilteredNames,
    noRows);
  data := transpose(data);
  matrixCol := size(data, 2);
  matrixRow := size(data, 1);
  noSteps := integer(data[matrixRow, 1] / timeSpan);

  dataDump := zeros(noSteps, matrixCol);
  i := 2; // iteration compares the previous row
  while i <= matrixRow loop
    if data[i - 1, 1] < timeSpan * timeStep then
      if data[i, 1] >= timeSpan * timeStep then
        dataDump[timeStep, :] := data[i, :];
        timeStep := timeStep + 1;
      end if;
    end if;
    i := i + 1;
  end while;

  DymolaCommands.MatrixIO.writeMatrix(
    outputFilename,
    "Environment",
    dataDump[:, 1:2]);

  DymolaCommands.MatrixIO.writeMatrix(
    outputFilename,
    "Zones",
    dataDump[:, 3:zoneCount], 
    true);

  DymolaCommands.MatrixIO.writeMatrix(
    outputFilename,
    "Walls",
    dataDump[:, zoneCount+1:varCount],
    true); // set true to append sheet to the same file

end MatrixConverter;
