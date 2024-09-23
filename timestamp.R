# Load necessary library
library(lubridate)

# Define log file path
log_file <- "log.csv"

# Function to log the current time
log_time <- function(log_file) {
  # Get the current time
  current_time <- Sys.time()
  
  # Write the time to the log file (append mode)
  write.table(data.frame(Time = current_time), 
              file = log_file, 
              sep = ",", 
              row.names = FALSE, 
              col.names = !file.exists(log_file), # add header only if the file doesn't exist
              append = TRUE)
}

cat("COUCOU")

# Infinite loop to log the time every minute
while (TRUE) {
  log_time(log_file)
  
  # Wait for 60 seconds (1 minute) before logging the next entry
  Sys.sleep(60)
}