#original code 
sp_P29182<-"clncpevycrrfqkc"
sp_P25668<-"cfitpd---itskdc"
sp_P01391<-"cfitpd---itskdc"
sp_P01382<-"cyktpip--itsetc"
hydropathy_scale<- c("i"=0.05,"v"=0.10,"l"=0.15,"f"=0.20,"c"=0.25,"m"=0.30,"a"=0.35,"g"=0.40,"t"=0.45,"s"=0.50,"w"=0.55,"y"=0.60,"p"=0.65,"h"=0.70,"q"=0.75,"n"=0.80,"e"=0.85,"d"=0.90,"k"=0.95,"r"=1.0)
hydropathy1 <- sapply(strsplit(sp_P01391, NULL)[[1]], function(aa) hydropathy_scale[aa])
hydropathy2 <- sapply(strsplit(sp_P01382, NULL)[[1]], function(aa) hydropathy_scale[aa])
calculate_rmsd <- function(values1, values2) {
  valid_indices <- !is.na(values1) & !is.na(values2)
  sqrt(mean((values1[valid_indices] - values2[valid_indices])^2))
}

# Calculate initial RMSD
initial_rmsd <- calculate_rmsd(hydropathy1, hydropathy2)
cat("Initial RMSD:", initial_rmsd, "\n")
optimize_rmsd_with_gaps <- function(hydro1, hydro2) {
  min_rmsd <- calculate_rmsd(hydro1, hydro2)
  best_hydro1 <- hydro1
  
  for (i in 1:(length(hydro1) - 1)) {
    # Skip swap if either position has a gap (NA)
    if (is.na(hydro1[i]) || is.na(hydro1[i + 1])) next
    
    # Swap adjacent residues
    swapped <- hydro1
    temp <- swapped[i]
    swapped[i] <- swapped[i + 1]
    swapped[i + 1] <- temp
    
    # Calculate new RMSD
    new_rmsd <- calculate_rmsd(swapped, hydro2)
    
    # Check if the new RMSD is lower
    if (new_rmsd < min_rmsd) {
      min_rmsd <- new_rmsd
      best_hydro1 <- swapped
    }
  }
  
  list(opt_hydro1 = best_hydro1, min_rmsd = min_rmsd)
}

# Run optimization
result <- optimize_rmsd_with_gaps(hydropathy1, hydropathy2)
cat("Optimized RMSD:", result$min_rmsd, "\n")
hydropathy_to_aa <- names(hydropathy_scale)[match(result$opt_hydro1, hydropathy_scale)]
final_sequence1 <- paste(hydropathy_to_aa, collapse = "")

cat("Optimized Sequence 1 Alignment:", final_sequence1, "\n")
cat("Sequence 2 Alignment:", sp_P01382, "\n")

--------------------------------------------------------------------------------
  #code update 
  # Kyte-Doolittle hydropathy scale
hydropathy_scale <- c(
    "i" = 0.05, "v" = 0.10, "l" = 0.15, "f" = 0.20, "c" = 0.25, "m" = 0.30, "a" = 0.35, "g" = 0.40, 
    "t" = 0.45, "s" = 0.50, "w" = 0.55, "y" = 0.60, "p" = 0.65, "h" = 0.70, "q" = 0.75, "n" = 0.80, 
    "e" = 0.85, "d" = 0.90, "k" = 0.95, "r" = 1.0
)

# Function to calculate RMSD
calculate_rmsd <- function(seq1, seq2) {
    # Remove gaps from both sequences
    seq1 <- gsub("-", "", seq1)
    seq2 <- gsub("-", "", seq2)
    
    # Ensure both sequences are of the same length after removing gaps
    len <- min(nchar(seq1), nchar(seq2))
    seq1 <- substr(seq1, 1, len)
    seq2 <- substr(seq2, 1, len)
    
    # Calculate hydropathy values for both sequences
    hydropathy1 <- sapply(strsplit(seq1, NULL)[[1]], function(x) hydropathy_scale[x])
    hydropathy2 <- sapply(strsplit(seq2, NULL)[[1]], function(x) hydropathy_scale[x])
    
    # Calculate RMSD
    rmsd <- sqrt(mean((hydropathy1 - hydropathy2)^2))
    return(rmsd)
}

# Function to optimize alignment by shifting gaps
optimize_alignment <- function(seq1, seq2, max_shifts = 10) {
    best_rmsd <- Inf
    best_alignment <- list(seq1, seq2)
    
    # Try all possible alignments by shifting gaps in both sequences
    for (shift1 in -max_shifts:max_shifts) {
        for (shift2 in -max_shifts:max_shifts) {
            # Shift seq1 and seq2 by adding '-' in front or at the end
            shifted_seq1 <- paste(rep("-", abs(shift1)), collapse = "")
            shifted_seq2 <- paste(rep("-", abs(shift2)), collapse = "")
            
            if (shift1 < 0) {
                shifted_seq1 <- paste0(shifted_seq1, seq1)
            } else {
                shifted_seq1 <- paste0(seq1, shifted_seq1)
            }
            
            if (shift2 < 0) {
                shifted_seq2 <- paste0(shifted_seq2, seq2)
            } else {
                shifted_seq2 <- paste0(seq2, shifted_seq2)
            }
            
            # Calculate RMSD for this alignment
            current_rmsd <- calculate_rmsd(shifted_seq1, shifted_seq2)
            
            # If RMSD is better, update the best alignment
            if (current_rmsd < best_rmsd) {
                best_rmsd <- current_rmsd
                best_alignment <- list(shifted_seq1, shifted_seq2)
            }
        }
    }
    
    # Return the best alignment
    return(best_alignment)
}

# Sample sequences (with gaps)
seq1 <- "clncpevycrrfqkc"
seq2 <- "cfitpd---itskdc"

# Optimize the alignment
optimized_alignment <- optimize_alignment(seq1, seq2)

# Print the optimized alignment and the RMSD
cat("Optimized Alignment:\n")
cat(paste(optimized_alignment[[1]], "\n"))
cat(paste(optimized_alignment[[2]], "\n"))
cat("Minimum RMSD:", calculate_rmsd(optimized_alignment[[1]], optimized_alignment[[2]]), "\n")
--------------------------------------------------------------------------------
  #code update - calculates the RMSD between the overall hydropathy 
  calculate_hydropathy <- function(seq) {
    # Map the amino acids to hydropathy values, ignore gaps (denoted by '-')
    hydropathy_values <- sapply(strsplit(seq, NULL)[[1]], function(aa) {
      if (aa == "-") {
        return(NA)  # Treat gaps as NA
      }
      return(hydropathy_scale[aa])
    })
    return(hydropathy_values)
  }

# Function to calculate RMSD between two sequences
calculate_rmsd <- function(seq1, seq2) {
  # Calculate hydropathy values for both sequences
  hydropathy1 <- calculate_hydropathy(seq1)
  hydropathy2 <- calculate_hydropathy(seq2)
  
  # Remove positions where either sequence has a gap (NA values)
  valid_indices <- which(!is.na(hydropathy1) & !is.na(hydropathy2))
  
  if (length(valid_indices) == 0) {
    return(NA)  # Return NA if no valid positions
  }
  
  # Compute the RMSD
  rmsd <- sqrt(mean((hydropathy1[valid_indices] - hydropathy2[valid_indices])^2))
  return(rmsd)
}

# Example usage:
seq1 <- "clncpevycrrfqkc"
seq2 <- "cfitp---ditskdc" 
rmsd_value <- calculate_rmsd(seq1, seq2)
print(paste("RMSD:", rmsd_value))
--------------------------------------------------------------------------------
#code update - calculates the RMSD in a way that it takes each pair of 
#residue in the alignment 
#handles residue alignment with a gap 
hydropathy_scale <- c("i"=0.05,"v"=0.10,"l"=0.15,"f"=0.20,"c"=0.25,
                      "m"=0.30,"a"=0.35,"g"=0.40,"t"=0.45,"s"=0.50,
                      "w"=0.55,"y"=0.60,"p"=0.65,"h"=0.70,"q"=0.75
                      ,"n"=0.80,"e"=0.85,"d"=0.90,"k"=0.95,"r"=1.0)

# Function to calculate hydropathy values for a sequence
calculate_hydropathy <- function(seq) {
  # Map the amino acids to hydropathy values, treat gaps ('-') as zero hydropathy
  hydropathy_values <- sapply(strsplit(seq, NULL)[[1]], function(aa) {
    if (aa == "-") {
      return(0)  # Treat gaps as hydropathy of 0
    }
    return(hydropathy_scale[aa])
  })
  return(hydropathy_values)
}

# Function to calculate RMSD between two aligned sequences
calculate_rmsd <- function(seq1, seq2) {
  # Calculate hydropathy values for both sequences
  hydropathy1 <- calculate_hydropathy(seq1)
  hydropathy2 <- calculate_hydropathy(seq2)
  
  # Remove positions where both sequences have gaps (0 values)
  valid_indices <- which(!(hydropathy1 == 0 & hydropathy2 == 0))
  
  if (length(valid_indices) == 0) {
    return(NA)  # Return NA if no valid positions
  }
  
  # Calculate RMSD for the corresponding residues in the alignment
  rmsd_sum <- sum((hydropathy1[valid_indices] - hydropathy2[valid_indices])^2)
  rmsd <- sqrt(rmsd_sum / length(valid_indices))  # RMSD formula
  return(rmsd)
}

# Example usage:
seq1 <- "cfkrf--nril---gkrydlgc"
seq2 <- "cytktwcdafcsirgkrvdlgc" 
rmsd_value <- calculate_rmsd(seq1, seq2)
print(paste("RMSD:", rmsd_value))

  
  
