library(tidyverse)

#######################################################################
## SQL for extracting manta tow data from oracle database            ##
## NOTE, the java applet 'dbExport.jar' is not part of this git      ##
## repository as it can only be successfully used by AIMS staff with ##
## certain granted privileges.                                       ##
## This script generates a flat file database extract (manta.csv).   ##
## Whilst not provided in this repository, it is available from AIMS ##
## and is linked in the supplementary materials.                     ##
#######################################################################
writeLines("
SELECT V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME,V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_LONG, V_RM_SAMPLE.REEF_LAT, 
  V_RM_SAMPLE.REPORT_YEAR, RM_MANTA.TOW_SEQ_NO, RM_MANTA.LIVE_CORAL, V_RM_SAMPLE.SAMPLE_CLASS
FROM RM_MANTA INNER JOIN V_RM_SAMPLE ON RM_MANTA.SAMPLE_ID = V_RM_SAMPLE.SAMPLE_ID
WHERE (((V_RM_SAMPLE.SAMPLE_CLASS) In ('K','C','G','Z') Or (V_RM_SAMPLE.SAMPLE_CLASS) Is Null))
ORDER BY V_RM_SAMPLE.P_CODE, V_RM_SAMPLE.A_SECTOR, V_RM_SAMPLE.SHELF, V_RM_SAMPLE.REEF_NAME, V_RM_SAMPLE.REEF_ID, V_RM_SAMPLE.REEF_LAT, V_RM_SAMPLE.REPORT_YEAR, 
RM_MANTA.TOW_SEQ_NO","data/primary/manta.sql")

if (goto_database_manta) system("java -jar scripts/dbExport.jar data/primary/manta.sql data/primary/manta.csv reef reefmon") 

manta <- read.csv('data/primary/manta.csv',strip.white=TRUE)  

save(manta, file='data/primary/manta.RData')

