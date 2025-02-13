import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# Load and process data
data = pd.read_csv("complexity_data.csv")
data = data.sort_values(by="d")


# Create plot
plt.figure(figsize=(10, 6))

# Calculate separate max values for x and y axes
x_max = max(data['d'])
y_max = max(max(data['C_complexity']), max(data['Rust_complexity']))  # Add 10% margin

# Create triangular regions using x_max
x = np.linspace(0, y_max, 2)
y = x

# Fill triangular regions
plt.fill_between(x, 0, y, color='green', alpha=0.1, label='Sub-linear')
plt.fill_between(x, y, y_max, color='red', alpha=0.1, label='Super-linear')

# Interpolate from 0 to first point
plt.plot([0, data['d'].iloc[0]], [0, data['C_complexity'].iloc[0]], 
         ':', color="#00aaff", linewidth=2)
plt.plot([0, data['d'].iloc[0]], [0, data['Rust_complexity'].iloc[0]], 
         ':', color="#ff1100", linewidth=2)

# Plot performance data
plt.plot(data['d'], data['C_complexity'], 'o-', color="#00aaff", 
         label="C Implementation", markersize=0, linewidth=2)
plt.plot(data['d'], data['Rust_complexity'], 'x-', color="#ff1100", 
         label="Rust Implementation", markersize=0, linewidth=2)

# Formatting
plt.grid(True, linestyle='--', alpha=0.5)
plt.title("Algorithm Performance: C vs Rust", fontsize=14)
plt.xlabel("String length $n$", fontsize=12)
plt.ylabel("Complexity $O(f(n))$", fontsize=12)
plt.legend(fontsize=10)

# Set different axis limits
plt.xlim(0, 300)
plt.ylim(0, y_max)

# Scientific notation
#plt.gca().xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
#plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
#plt.ticklabel_format(style="sci", axis="both", scilimits=(0,0))

plt.tight_layout()
plt.savefig("complexity_comparison.png", dpi=300)
plt.show()