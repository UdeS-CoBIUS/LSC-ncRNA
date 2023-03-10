import matplotlib.pyplot as plt

# Read in memory usage data from log file
with open('memory.log', 'r') as f:
    lines = f.readlines()

# Extract memory usage values and timestamps
memory_usage = []
timestamps = []
for line in lines:
    if 'Memory usage' in line:
        memory_usage.append(float(line.split(': ')[-1].strip().replace('GB', '')))
        timestamps.append(line.split(' ')[0] + ' ' + line.split(' ')[1])

# Create line plot
plt.plot(timestamps, memory_usage)
plt.title('Memory Usage Over Time')
plt.xlabel('Time')
plt.ylabel('Memory Usage (GB)')
plt.xticks(rotation=45)
plt.show()
