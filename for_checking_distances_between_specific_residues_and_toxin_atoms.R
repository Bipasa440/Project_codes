library(bio3d)

# Read the PDB file
pdb <- read.pdb("C:/Users/Bipasa/Downloads/APU_SEM_6/Honours_work/for_docking/1yi5_1.pdb")

# Select all atoms of the 89th residue in chain A (tyrosine)
chain_a_89 <- atom.select(pdb, chain = "A", resno = 89)

# Ensure the selected residue is tyrosine
if (unique(pdb$atom[chain_a_89$atom, "resid"]) != "TYR") {
  stop("Residue 89 in chain A is not tyrosine.")
}

# Extract atom coordinates of chain A residue 89
atoms_a_89 <- pdb$atom[chain_a_89$atom, ]

# Select all atoms in chain C
chain_c <- atom.select(pdb, chain = "C")
atoms_c <- pdb$atom[chain_c$atom, ]

# Function to calculate pairwise distances without memory overhead
calculate_distances_efficient <- function(atoms1, atoms2) {
  result <- data.frame()  # Initialize an empty data frame to store results
  
  # Loop through each atom in atoms1
  for (i in seq_len(nrow(atoms1))) {
    coord1 <- as.numeric(atoms1[i, c("x", "y", "z")])  # Ensure coord1 is numeric
    
    # Calculate distances to all atoms in atoms2
    distances <- sqrt(rowSums((t(t(atoms2[, c("x", "y", "z")]) - coord1))^2))
    
    # Create a data frame for the current atom
    temp <- data.frame(
      chain1 = atoms1[i, "chain"],
      res1 = atoms1[i, "resno"],
      resid1 = atoms1[i, "resid"],
      atom1 = atoms1[i, "elety"],
      chain2 = atoms2$chain,
      res2 = atoms2$resno,
      resid2 = atoms2$resid,
      atom2 = atoms2$elety,
      distance = round(distances, 2)  # Round distances to 2 decimal places
    )
    
    # Append results to the main data frame
    result <- rbind(result, temp)
  }
  
  return(result)
}

# Calculate distances
cat("Calculating distances...\n")
start_time <- Sys.time()
distances <- calculate_distances_efficient(atoms_a_89, atoms_c)
end_time <- Sys.time()
cat("Distance calculations completed in", difftime(end_time, start_time, units = "secs"), "seconds.\n")

# Display the distances as a data frame
cat("Calculated distances:\n")
print(distances)
#filteration
filtered_distances <- distances[distances$distance < 6, ]
#this confirms what I found - so is chimera wrong (but how and why)
