import numpy as np
import matplotlib.pyplot as plt
import scipy.io

# Parameters
total_length = 3.0  # meters
transition_depth = 1.0
dz_min = 0.0022
dz_max = 0.1818
uniform_dz = dz_max
total_nodes = 36
total_intervals = total_nodes - 1

# Calculate uniform intervals
uniform_length = total_length - transition_depth
uniform_intervals = int(np.round(uniform_length / uniform_dz))

# Logarithmic intervals
log_intervals = total_intervals - uniform_intervals

# Logarithmic spacing and nodes
log_spacings = np.geomspace(dz_min, dz_max, num=log_intervals)



log_nodes = np.concatenate(([0], np.cumsum(log_spacings)))
log_nodes *= (transition_depth / log_nodes[-1])  # Scale to 1 m

# Uniform spacing and nodes
uniform_nodes = np.arange(log_nodes[-1] + uniform_dz, total_length + 1e-6, uniform_dz)

# All nodes
all_nodes = np.concatenate((log_nodes, uniform_nodes))

# Element sizes (depth spacing between nodes)
dz_all = np.diff(all_nodes)

# Midpoints of elements for plotting
midpoints = (all_nodes[:-1] + all_nodes[1:]) / 2
node_indices = np.arange(1, total_nodes)


# Recalculate node centers and corresponding node numbers
node_centers = (all_nodes[:-1] + all_nodes[1:]) / 2
node_numbers = np.arange(1, len(node_centers) + 1)


# Plot
# Recalculate all required variables for plotting
# dz_all = np.diff(all_nodes)
# midpoints = (all_nodes[:-1] + all_nodes[1:]) / 2
# node_numbers = np.arange(1, len(midpoints) + 1)
# node_centers = midpoints

# Create 1 by 2 subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6), sharey=True)

# First plot: Depth vs Node Size
ax1.plot(dz_all, midpoints, marker='o', linestyle='-')
ax1.invert_yaxis()
ax1.set_xlabel("Node size (m)")
ax1.set_ylabel("Depth (m)")
ax1.set_title("Depth vs Node Size")
ax1.grid(True)

# Second plot: Depth vs Centered Node Number
ax2.plot(node_numbers, node_centers, marker='o', linestyle='-')
ax2.invert_yaxis()
ax2.set_xlabel("Node Number (centered in each block)")
ax2.set_title("Depth vs Centered Node Number")
ax2.set_xlim(0, 40)
ax2.set_ylim(0, 3)
ax2.grid(True)

# Layout adjustment
plt.tight_layout()
plt.savefig("grid_layout.png",dpi = 200)


# Save variables to a .mat file
mat_data = {
    'dz_all': dz_all,
    'z_req': midpoints
}

scipy.io.savemat('grid_spacing.mat', mat_data)
















