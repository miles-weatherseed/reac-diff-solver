import matplotlib.pyplot as plt
import matplotlib.animation as animation

def animate(u, total_time):
    """
    Creates an animation of u
    :param u: frames to show
    :param total_time: total time over which the animation should run
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
    plt.show()
