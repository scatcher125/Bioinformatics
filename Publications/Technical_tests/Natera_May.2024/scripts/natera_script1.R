# R script

# load libraries
library(openxlsx)
library(dplyr)

# 1. load data
markers_df <- openxlsx::read.xlsx("~/Bioinformatics/Interviews/Natera_data.xlsx",
                                  sheet="Markers",
                                  detectDates = TRUE)
patients_df <- openxlsx::read.xlsx("~/Bioinformatics/Interviews/Natera_data.xlsx",
                                  sheet="Patients",
                                  detectDates = TRUE)
colnames(patients_df)[2] <- "patient_ID"

# 2. From the table Markers, removes all records with the "Blood draw date" before 2019
markers_df <- markers_df %>%
  dplyr::filter(Blood.draw.date < Date(2020-01-01))



# 3. Creates a third table that is a merge of the first two (i.e. contains all the columns from both tables), by Patient ID.
merged_df <- merge(markers_df,
                   patients_df,
                   by="patient_ID")


# 4. To the Patients table, adds a column showing how many CEA assays were taken for 
# each patient (for example it will be 2 for EN6K, 6 for FKB8, etc). Save the table as a tab-separated file.
write.table(patients_df,
            file="~/Bioinformatics/Interviews/patient_df.txt",
            sep="\t")

# 5. In the Markers table, adds a column with relative date, showing the number of days from surgery to this blood draw (for example, a sample taken on 03/29/2020 for patient "CHIQ" will have a relative date of 231 days because this patient had surgery on 08/10/2020).

# 6. For each patient, draws a plot showing how the levels of those markers changed over time
