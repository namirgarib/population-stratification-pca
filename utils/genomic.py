import matplotlib.pyplot as plt
import numpy as np
import random

# Configuration
SEGMENTS = 500  # Total number of color chunks per track (must be even)
CHUNK_SIZE = 20  # Number of pixels per color chunk
SEPARATOR_HEIGHT = 1  # Blank rows between tracks
COLORS = {
    'black': [0, 0, 0],
    'red': [255, 0, 0],
    'blue': [0, 0, 255],
    'green': [0, 255, 0]
}

def create_track(base_colors, new_color):
    """Create a new track with 50% previous colors and 50% new color"""
    n = len(base_colors)
    preserved = random.sample(base_colors, n // 2)
    new = [new_color] * (n // 2)
    combined = preserved + new
    random.shuffle(combined)
    return combined

# Generate tracks
track1 = ['black'] * SEGMENTS

track2 = ['black'] * (SEGMENTS//2) + ['red'] * (SEGMENTS//2)
random.shuffle(track2)

track3 = create_track(track2, 'green')
track4 = create_track(track3, 'blue')

# Build image with separators
image = []
tracks = [track1, track2, track3, track4]
track_indices = []

for i, track in enumerate(tracks):
    # Expand each color chunk to multiple pixels
    expanded_track = []
    for color in track:
        expanded_track.extend([COLORS[color]] * CHUNK_SIZE)
    
    # Add track to image
    image.append(expanded_track)
    track_indices.append(len(image)-1)  # Track position
    
    # Add separator if not last track
    if i < len(tracks) - 1:
        separator = [[255, 255, 255]] * len(expanded_track)  # White
        for _ in range(SEPARATOR_HEIGHT):
            image.append(separator)

# Convert to numpy array
image_array = np.array(image)

# Plotting
plt.figure(figsize=(10, 6))
plt.imshow(image_array, aspect='auto', interpolation='nearest')
plt.yticks(track_indices, ['GEN0', 'GEN1', 'GEN2', 'GEN3'])
plt.xticks([])
plt.tight_layout()
plt.show()