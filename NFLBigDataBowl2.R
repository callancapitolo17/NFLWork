# Specify the path to your zip file
zip_path <- "/Users/callancapitolo/NFLWork/nfl-big-data-bowl-2025.zip"

# List files in the zip archive
file_list <- unzip(zip_path, list = TRUE)$Name

# Filter to include only CSV files
csv_files <- file_list[grepl("\\.csv$", file_list)]

# Read each CSV file and store them in a list
csv_data_list <- lapply(csv_files, function(file) {
  # Read each CSV file from the zip
  read.csv(unz(zip_path, file))
})

# Optionally, name each element in the list after its corresponding CSV file
names(csv_data_list) <- csv_files

# Now csv_data_list contains a separate data frame for each CSV file in the zip
csv_data_list$games.csv
