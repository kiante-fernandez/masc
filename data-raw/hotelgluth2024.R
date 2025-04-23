## code to prepare `hotelgluth2024` dataset goes here
hotelgluth2024 <- merge(merge(merge(
  read.csv("~/Downloads/masc_cr_KF/data/attValues_Hotel_all_subjects_JD.csv"),
  read.csv("~/Downloads/extracted_hotel_data.csv"),
  all.x = TRUE
),
transform(
  t(apply(read.csv("~/Downloads/masc_cr_KF/data/attWeights_Hotel.csv"), 1, function(x) x / sum(x))),
  subject = read.csv("~/Downloads/masc_cr_KF/data/modelParams_Hotel_all_subjects.csv")$subject
),
all.x = TRUE
),
read.csv("~/Downloads/masc_cr_KF/data/modelParams_Hotel_all_subjects.csv"),
all.x = TRUE
)

# Add dataset column
hotelgluth2024$dataset <- "hotel"

# Reorder columns
hotelgluth2024 <- hotelgluth2024[,c("dataset","subject", "trial", "opt1_att1", "opt1_att2", "opt1_att3",
                                    "opt2_att1", "opt2_att2", "opt2_att3", "difficulty", "choice",
                                    "rt", "att_w1", "att_w2", "att_w3", "sigma", "alpha", "delta"
)]

# Sort by subject and trial
hotelgluth2024 <- hotelgluth2024[order(hotelgluth2024$subject, hotelgluth2024$trial), ]

usethis::use_data(hotelgluth2024, compress = "xz", overwrite = TRUE)
