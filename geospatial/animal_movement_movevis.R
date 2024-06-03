# Install and load necessary packages
install.packages("move")
library(move)

# Log in to Movebank
loginStored <- movebankLogin(username="yourusername", password="yourpass")

# Search for Movebank studies related to migration
searchMovebankStudies(x="migration", login=loginStored)

# Get Movebank ID for a specific study
getMovebankID("Ocelots on Barro Colorado Island, Panama", login=loginStored)

# Get data for a specific Movebank study
getMovebankStudy(study="Ocelots on Barro Colorado Island, Panama", login=loginStored)

# Get Movebank data for a specific study
swiss <- getMovebankData(study="Wheatear migration behaviour of an alpine breeding population (Switzerland)", login=loginStored)

# Install and load moveVis package
devtools::install_github("16EAGLE/moveVis")
library(moveVis)

# Align move_data to a uniform time scale
m <- align_move(swiss, res = 4, unit = "mins")

# Create spatial frames with an OpenStreetMap watercolour map
frames <- frames_spatial(m, path_colours = NA,
                         map_service = "esri",map_type = NULL, alpha = 0.5) %>% 
  add_labels(x = "Longitude", y = "Latitude") %>% # add some customizations, such as axis labels
  add_northarrow() %>% 
  add_scalebar() %>% 
  add_timestamps(type = "label") %>% 
  add_progress()

# Preview one of the frames, e.g. the 100th frame
frames[[100]]

# Animate frames and save as a GIF
animate_frames(frames, out_file = "youranimationname.gif")
