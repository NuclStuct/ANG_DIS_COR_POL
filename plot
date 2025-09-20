# Set terminal type and output file
set terminal pngcairo size 800, 600 enhanced font "Arial,10"
set output "q4.png"

# Add a title and labels
set title "Heatmap Example"
set xlabel "teta"
set ylabel "energy"

# Define color palette
set palette rgbformulae 33,13,10  # Or use "set palette defined" for custom colors
# Example custom palette: set palette defined (0 "blue", 0.5 "white", 1 "red")

# Set up color bar
set cblabel "Value"
set cbrange [0:2]  # Adjust based on your data range

# Configure 2D map view
set view map
unset key  # Disable legend
set grid   # Add grid

# Plot the data as an image
splot 'q4.dat' using 1:2:3 with image

