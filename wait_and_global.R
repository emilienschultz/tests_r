# Number of minutes to wait
minutes_to_wait <- 90

# Loop through each minute
for (i in 1:minutes_to_wait) {
  # Print the current time
  print(paste("Time:", Sys.time()))
  
  # Sleep for 1 minute (60 seconds)
  Sys.sleep(60)
}

# After the loop, create a new variable in the global environment
assign("script_ended", TRUE, envir = .GlobalEnv)

# Optionally, print a message indicating the script has ended
print("The script has ended and 'script_ended' has been created.")