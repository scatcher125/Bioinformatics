## TNYA0044 Format meta data
## Reva Shenwai
## 15-Jun-2024
## ---------------------------------------------------

# Libraries
library(dplyr)

# load data
meta_data <- read.csv(file="~/Bioinformatics/Interviews/Seminar/data/TNYA0044/TNYA0044_metadata_orig.csv")

# format metadata
meta_data$Sample_ID <- sapply(meta_data$Full_sample_ID, 
                              function(samp_id) {
                                vec <- unlist(strsplit(samp_id,"_"))
                                to_ret <- paste0("X",
                                       vec[length(vec)-1],
                                       "_",
                                       vec[length(vec)]
                                       )
                                return(to_ret)
                                })

meta_data$Group <- sapply(meta_data$Full_sample_ID, 
                          function(samp_id) {
                            vec <- unlist(strsplit(samp_id,"_"))
                            if (length(vec)==4) {
                              to_ret <- paste0(unlist(strsplit(samp_id,"_"))[1],
                                               "_",
                                               unlist(strsplit(samp_id,"_"))[2]
                                               )
                            } else {
                              to_ret <- paste0(unlist(strsplit(samp_id,"_"))[1],
                                               "_",
                                               unlist(strsplit(samp_id,"_"))[2],
                                               "_",
                                               unlist(strsplit(samp_id,"_"))[3]
                              )
                            }
                            })

# save file
write.csv(meta_data,
          file="~/Bioinformatics/Interviews/Seminar/data/TNYA0044/TNYA0044_metadata.csv")


