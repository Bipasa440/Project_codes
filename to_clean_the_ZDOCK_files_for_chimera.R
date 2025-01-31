# Define a function to truncate each line after the 54th character
clean_pdb_line <- function(line) {
  # Keep only the first 54 characters of the line
  cleaned_line <- substring(line, 1, 54)
  return(cleaned_line)
}

# Read the PDB file
input_file <- "C:/Users/Bipasa/Downloads/APU_SEM_6/Honours_work/for_docking/complex.1.pdb"
output_file <- "C:/Users/Bipasa/Downloads/APU_SEM_6/Honours_work/for_docking/complex.1_cleaned.pdb"
pdb_lines <- readLines(input_file)

# Apply the cleaning function to each line
cleaned_lines <- sapply(pdb_lines, clean_pdb_line)

# Write the cleaned lines to a new file
writeLines(cleaned_lines, output_file)

cat("The cleaned PDB file has been saved to", output_file, "\n")
