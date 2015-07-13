
#station= "Fichtelberg"
station= "KahlerAsten"

mode= "std"
#mode= "mcs"

mcs_nSamples= 200

dir_base= "/home/dkneis/progress/echse_apps/hypsoRR/processes/snow/modeltest"

dir_meteo= paste(dir_base,"/input/DWD/use/",station,sep="")
file_obs= paste(dir_base,"/input/DWD/raw/",station,"/snowwe.txt",sep="")
col_obs= "snowwe"
fact_obs= 0.001

col_time= "date"

file_timesGOF= paste(dir_base,"/input/timesGOF.txt",sep="") 

file_sensors= paste(dir_base,"/input/sensors.txt",sep="")
file_initials= paste(dir_base,"/input/initials.txt",sep="")

ofile_sim= paste(dir_base,"/output/",station,"/snow_sim.txt",sep="")
ofile_dbg= paste(dir_base,"/output/",station,"/snow_dbg.txt",sep="")
ofile_mcs_runs= paste(dir_base,"/output/",station,"/snow_mcs_runs.txt",sep="")
ofile_mcs_best= paste(dir_base,"/output/",station,"/snow_mcs_best.txt",sep="")
ofile_pdf= paste(dir_base,"/output/",station,"/snow_sim.pdf",sep="")
replace= TRUE

# Table of model parameters
file_params= paste(dir_base,"/input/paramRanges.txt",sep="")
#file_params= ofile_mcs_best

# Number of sub-time steps per time step
n_dtsub= 12

