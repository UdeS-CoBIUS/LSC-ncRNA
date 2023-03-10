import psutil
import time
import logging

# Configure logging
logging.basicConfig(filename='memory.log', level=logging.INFO)

# Define function to get memory usage and log it
def log_memory_usage():
    # Get memory usage in GB
    memory_usage = psutil.virtual_memory().used / 1024 / 1024 / 1024
    # Get current time
    current_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    # Log memory usage and timestamp to file
    logging.info(f"{current_time} Memory usage: {memory_usage:.2f} GB")

# Schedule task to run every 2.5 minutes
while True:
    log_memory_usage()
    time.sleep(150) # sleep for 2.5 minutes
