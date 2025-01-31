#loading packages and libraries
library(bio3d)
library(parallel)

# Read PDB file
pdb <- read.pdb("C:/Users/Bipasa/Downloads/APU_SEM_6/Honours_work/for_docking/1yi5_1.pdb")
cat("Extracting chains...\n") #debugging statement
chains <- list(
  A = atom.select(pdb, chain = "A", model = 1),
  B = atom.select(pdb, chain = "B", model = 1),
  C = atom.select(pdb, chain = "C", model = 2) # Model 2 for chain C
)
cat("Chains extracted.\n")

# Extract atoms for each chain using lapply and atom.select()
atoms <- lapply(chains, function(chain) pdb$atom[chain$atom, ])

# Function to calculate distances and filter them
calculate_distances <- function(atoms1, atoms2, threshold = 6) {
  # Check if either input is empty, if not then proceed (debugging statement)
  if (nrow(atoms1) == 0 || nrow(atoms2) == 0) {
    return(data.frame())  # Return empty data frame if no atoms
  }
  #extract coordinates
  coords1 <- atoms1[, c("x", "y", "z")]
  coords2 <- atoms2[, c("x", "y", "z")]
  
  # Vectorized calculation of pairwise distances
  distances <- as.matrix(dist(rbind(coords1, coords2)))
  #subsetting such that distances within the chain is ignored
  distances <- distances[1:nrow(coords1), (nrow(coords1) + 1):nrow(distances)]
  
  # Filter by threshold (distance <= threshold)
  filtered <- which(distances <= threshold, arr.ind = TRUE)
  
  # If no distances meet the threshold, return an empty data frame
  if (length(filtered) == 0) {
    return(data.frame())
  }
  
  # Store results in a data frame
  result <- data.frame(
    chain1 = atoms1[filtered[, 1], "chain"],
    res1 = atoms1[filtered[, 1], "resno"],
    resid1 = atoms1[filtered[, 1], "resid"],
    atom1 = atoms1[filtered[, 1], "elety"],
    chain2 = atoms2[filtered[, 2], "chain"],
    res2 = atoms2[filtered[, 2], "resno"],
    resid2 = atoms2[filtered[, 2], "resid"],
    atom2 = atoms2[filtered[, 2], "elety"],
    distance = signif(distances[filtered], 2)
  )
  
  return(result)
}

# Define chain pairs for comparison
chain_pairs <- list(
  list(atoms$A, atoms$C),
  list(atoms$B, atoms$C)
)

# Perform parallel distance calculations
cat("Calculating distances in parallel...\n")
start_time <- Sys.time()

# Detect the number of available cores and create a cluster, leave one 
#core for other tasks
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

# Export required data and functions to the cluster, export to each core
clusterExport(cl, c("calculate_distances", "atoms"))
clusterEvalQ(cl, library(bio3d)) # Load necessary libraries on each worker

# Apply calculations in parallel
results_list <- parLapply(cl, chain_pairs, function(pair) {
  calculate_distances(pair[[1]], pair[[2]], threshold = 6)
})

# Stop the cluster after calculations
stopCluster(cl)

end_time <- Sys.time()
cat("Distance calculations completed in", difftime(end_time, start_time, units = "secs"), "seconds.\n")

# Combine results
results <- do.call(rbind, results_list)

# Display and save results
cat("Results:\n")
print(results)

#Filtering out the residue pairs 
residue_pairs <- unique(results[, c("chain1", "res1", "resid1", "chain2", "res2", "resid2")])

# Display the extracted residue pairs
cat("Filtered Residue Pairs:\n")
print(residue_pairs)

#save as csv (optional)
write.csv(results, "C:/Users/Bipasa/Downloads/APU_SEM_6/Honours_work/for_docking/optimised_code_results.csv", row.names = FALSE)