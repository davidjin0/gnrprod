library(tidyverse)
library(haven)
library(pracma)

data_311 <- read_dta('Replication_files/Tables_2_3/Colombia/311/data_col.dta')
industry_311 <- prodgnr(RGO ~ L + K + RI, share = si, id = id, time = year,
                        data = data_311)

data_321 <- read_dta('Replication_files/Tables_2_3/Colombia/321/data_col.dta')
industry_321 <- prodgnr(RGO ~ L + K + RI, share = si, id = id, time = year,
                        data = data_321)

data_322 <- read_dta('Replication_files/Tables_2_3/Colombia/322/data_col.dta')
industry_322 <- prodgnr(RGO ~ L + K + RI, share = si, id = id, time = year,
                        data = data_322)

data_331 <- read_dta('Replication_files/Tables_2_3/Colombia/331/data_col.dta')
industry_331 <- prodgnr(RGO ~ L + K + RI, share = si, id = id, time = year,
                        data = data_331)

data_381 <- read_dta('Replication_files/Tables_2_3/Colombia/381/data_col.dta')
industry_381 <- prodgnr(RGO ~ L + K + RI, share = si, id = id, time = year,
                        data = data_381)







