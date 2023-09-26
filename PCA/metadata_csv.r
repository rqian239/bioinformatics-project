library(stringr)

determine_subject_type <- function(subject_string) {
    if (str_detect(subject_string, "bipolar disorder")) {
        return("bipolar disorder")
    } else if (str_detect(subject_string, "major depression")) {
        return("depression")
    } else if (str_detect(subject_string, "schizophrenia")) {
        return("schizophrenia")
    } else {
        return("control")
    }
}

metadata <- read.csv(file = "./PCA/metadata_PCA.tsv", sep = "\t")

for (i in 1:nrow(metadata)) {
    metadata[i,]$refinebio_disease <- determine_subject_type(metadata[i,]$refinebio_subject)
}

write.csv(metadata, file = "metadata_PCA.csv")