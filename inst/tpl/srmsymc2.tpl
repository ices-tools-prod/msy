DATA_SECTION
 !! ad_comm::change_datafile_name("srmsymc.par");
  init_number ap
  init_number bp
 LOCAL_CALCS
  ofstream srmsymc_pe("srmsymc_pe.dat");
  srmsymc_pe<<ap<<" "<<bp<<endl;
  srmsymc_pe.close();
 END_CALCS

PARAMETER_SECTION
  objective_function_value f;

PROCEDURE_SECTION
