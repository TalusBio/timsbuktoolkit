# /// script
# requires-python = ">=3.12"
# dependencies = [
#    "matplotlib",
# ]
# ///

from matplotlib import pyplot as plt
import json

# Side by side plot
left = "centroided_peaks.json"
right = "input_peaks.json"  # Each is a list[(float, float, float)] where each tuple is (mz, ims, intensity)

left_data = json.load(open(left))
right_data = json.load(open(right))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True)
min_inten = min(p[2] for p in left_data)
ax1.scatter(
    [p[0] for p in left_data],
    [p[1] for p in left_data],
    s=[p[2] / 1000 for p in left_data],  # Scale down intensity for visibility
    # c=[0.5 if p[2] == min_inten else 1 for p in left_data],  # Color by intensity
    alpha=0.5,
)

min_inten = min(p[2] for p in right_data)

ax2.scatter(
    [p[0] for p in right_data],
    [p[1] for p in right_data],
    s=[p[2] / 1000 for p in right_data],  # Scale down intensity for visibility
    # c=[0.5 if p[2] == min_inten else 1 for p in right_data],  # Color by intensity
    alpha=0.5,
)

fig.tight_layout()
plt.show()
