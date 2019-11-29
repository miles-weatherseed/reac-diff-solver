import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter


def animate(u, total_time, filename = ""):
    """
    Creates an animation of u
    :param u: frames to show
    :type u: list of 2d numpy arrays
    :param total_time: total time over which the animation should run
    :type total_time: float
    :param filename: file to which to save the animation, if no filename is provided, the animation will not be saved
    :type filename: string
    :return:
    """
    fig = plt.figure()

    ims = []

    for frame in u:
        im = plt.imshow(frame, animated=True)
        ims.append([im])

    interval = (1000*total_time) / len(u) # number of miliseconds between each frame
    ani = animation.ArtistAnimation(fig, ims, interval=interval, blit=True,
                                    repeat_delay=1000)

    if filename != "":
        writer = PillowWriter(fps=15, bitrate=1800)
        ani.save(filename + ".gif", writer = writer)

    plt.show()
