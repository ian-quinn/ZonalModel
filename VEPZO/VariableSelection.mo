within VEPZO;
model VariableSelection
 annotation (__Dymola_selections={Selection(name="ZoneTemp", match={MatchVariable(name="*.medium.T")})});
end VariableSelection;
