# This script will format data into stylish tables
# From https://www.littlemissdata.com/blog/prettytables

# # Install packages (do this once)

# install.packages("data.table")
# install.packages("dplyr")
# install.packages("formattable")
# install.packages("tidyr")

#Load the libraries

library(data.table)
library(dplyr)
library(formattable)
library(tidyr)

#Set a few color variables to make our table more visually appealing

customGreen0 = "#DeF7E9"
customGreen = "#71CA97"
customRed = "#ff7f7f"

# Where the data is!
data_file_path <- "./results/GO_results.tsv"

# Read in the data
main_df <- fread(data_file_path, sep = "\t", data.table = FALSE, header = TRUE, stringsAsFactors = FALSE)

# View the table in RStudio, you can export it as an image
formattable(main_df)