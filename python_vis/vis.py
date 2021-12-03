import numpy as np
import matplotlib.pyplot as plt
import difflib
import matplotlib.animation as animation
import matplotlib.cm as cm


with open('rho_test2.txt', 'r') as f:
    dat = f.readlines()

fig = plt.figure()

# make_image = lambda img: [plt.imshow(np.array(list(map(float, img.split()))).reshape(256, 256).swapaxes(0, 1), cmap='viridis', animated=True)]
make_image = lambda img: [plt.imshow(np.array(list(map(float, img.split()))).reshape(256, 256), cmap='viridis', animated=True)]


frames = [make_image(img) for img in dat]

ani = animation.ArtistAnimation(fig, frames, interval=50, blit=True)

writer = animation.FFMpegWriter(
    fps=15, metadata=dict(artist='Me'), bitrate=1800)
ani.save("movie2.mp4", writer=writer)
