#' title: "Protein aggregate area analysis"

#' author: "Bipasa Show"

#' date: "27/10/24"

#Install BiocManager and EMImage before doing all these. 
#Also make sure your image is in your workig directory 

# Load and process image as before
img <- readImage("C:/Users/Bipasa/Downloads/APU SEM 5/Honours_work/
                 aamodis_stuff/1_cl.jpg")  # Replace with the actual path of your image


# Isolate the green channel (assuming RGB image)
green_channel <- channel(img, "green")

# Enhance contrast to better isolate green regions (optional, adjust as necessary)
green_channel <- normalize(green_channel)

# Thresholding - adjust threshold to capture more meaningful green areas
binary_img <- green_channel > 0.5  # Adjust threshold value as needed

# Optional: Morphological operations to clean up noise with a smaller square brush
binary_img <- opening(binary_img, makeBrush(2, shape = "box"))
binary_img <- closing(binary_img, makeBrush(2, shape = "box"))

# Label the objects (aggregates) in the binary image
label_img <- bwlabel(binary_img)

# Calculate the area of each aggregate
aggregate_areas <- computeFeatures.shape(label_img)[, "s.area"]

# Create a data frame of all aggregates without filtering by area
area_table <- data.frame(Aggregate_ID = seq_along(aggregate_areas),
                         Area_in_Pixels = aggregate_areas)

# Display the original image with labeled aggregates
display(img, method = "raster")
title("Protein Aggregates with IDs")

# Overlay labels for each aggregate
for (i in 1:nrow(area_table)) {
  # Find the coordinates for the aggregate's centroid
  aggregate_id <- area_table$Aggregate_ID[i]
  centroid <- computeFeatures.moment(label_img)[aggregate_id, c("m.cx", "m.cy")]
  
  # Add text label at the centroid
  text(centroid[1], centroid[2], labels = aggregate_id, col = "red", cex = 0.6)
}

# Print the area table
print(area_table)