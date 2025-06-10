import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

def AnimateGrid2d(bounds, update_method, threeColors=False, frames=100, interval=50, scale=1, extremis=None, cmap='winter', save_as=None):
    """
    Creates an animation of a 2D grid.

    Parameters:
    - bounds: tuple of (x, y) dimensions for the grid
    - update_method: object with methods `module_stepper`, `module_call`, and `module_title`
    - threeColors: bool, whether to use three colors (RGB) for the grid
    - frames: int, number of frames in the animation
    - interval: int, time interval between frames in milliseconds
    - scale: int, scaling factor for the grid
    - extremis: tuple of (vmin, vmax) for color scale limits
    - cmap: str, colormap for the grid
    - save_as: str, file name to save the animation (e.g., 'animation.mp4')
    """
    if threeColors:
        fig, ax = plt.subplots()

        nx = bounds[0]
        ny = bounds[1]

        def update(frame):
            update_method.module_stepper()
            ax.clear()

            # Generate RGB values for each cell
            rgb_grid = np.zeros((ny * scale, nx * scale, 3))

            for j in range(ny * scale):
                for i in range(nx * scale):
                    rgb_grid[j,i] = update_method(i % nx, j % ny)
            # Plot the grid
            ax.imshow(rgb_grid, origin='lower', extent=[0, nx * scale, 0, ny * scale], interpolation='nearest')
            ax.set_title(update_method.module_title())
    else:
        fig, (ax_top, ax_bottom) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [10, 1]})

        nx = bounds[0]
        ny = bounds[1]

        def update(frame):
            update_method.module_stepper()
            ax_top.clear()
            ax_bottom.clear()

            # Generate RGB values for each cell
            rgb_grid = np.zeros((ny * scale, nx * scale))  # Rows = y, Cols = x

            max_values = []
            for j in range(ny * scale):
                maxi = 0
                for i in range(nx * scale):
                    rgb_grid[j, i] = update_method(i % nx, j % ny)
                    if i == nx * scale // 2:
                        maxi = rgb_grid[j, i]

                max_values.append(maxi)
            # Plot the grid
            if extremis:
                ax_top.imshow(rgb_grid, origin='lower', extent=[0, nx * scale, 0, ny * scale], cmap=cmap, interpolation='nearest', 
                              vmin=extremis[0], vmax=extremis[1])
            else:
                ax_top.imshow(rgb_grid, origin='lower', extent=[0, nx * scale, 0, ny * scale], cmap=cmap, interpolation='nearest')
            ax_bottom.plot(range(ny * scale), max_values, color='blue')
            ax_bottom.set_xlim(0, ny * scale)
            ax_bottom.set_ylim(min(max_values) * 1.1, max(max_values) * 1.1)
            ax_top.set_title(update_method.module_title())

    ani = FuncAnimation(fig, update, frames=frames, interval=interval)

    if save_as:
        ani.save(save_as, writer='ffmpeg')
        print(f'Animation saved as {save_as}')
    else:
        plt.show()

# Example usage:
# AnimateGrid2d((10, 10), update_method, save_as='animation.mp4')

def AnimateScatter(bounds, stepper, positions, frames=100, interval=50, save_as=None,
                   get_color_field=None, cmap='viridis', color_limits=None, background=None):
    """
    Animate a scatter plot of points, optionally colored by a field value.

    Parameters:
    - bounds: (xmin, xmax, ymin, ymax)
    - stepper: object with Step() method
    - positions: object with .size() and indexable Vector-like values (.x, .y)
    - frames: number of animation frames
    - interval: milliseconds between frames
    - save_as: filename to save animation
    - get_color_field: callable i -> float value, used to color each point
    - cmap: matplotlib colormap name
    - color_limits: tuple (vmin, vmax) for color normalization
    - background: optional color string or RGB tuple (sets axes and figure background)
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.animation import FuncAnimation
    from matplotlib.colors import Normalize

    fig, ax = plt.subplots()
    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[2], bounds[3])

    if background:
        fig.patch.set_facecolor(background)
        ax.set_facecolor(background)

    scat = ax.scatter([], [], c=[], cmap=cmap)
    norm = Normalize(*color_limits) if color_limits else None

    def init():
        scat.set_offsets(np.empty((0, 2)))
        return scat,

    def update(frame):
        stepper.Step()

        points = [positions[i] for i in range(positions.size())]
        xy = np.array([[p.x, p.y] for p in points])
        scat.set_offsets(xy)

        if get_color_field:
            color_values = np.array([get_color_field(i) for i in range(positions.size())])
            scat.set_array(color_values)
            if color_limits:
                scat.set_clim(*color_limits)

        return scat,

    ani = FuncAnimation(fig, update, init_func=init, frames=frames, interval=interval, blit=True)

    if save_as:
        ani.save(save_as, writer='ffmpeg')
        print(f'Animation saved as {save_as}')
    else:
        plt.show()



class AnimationUpdateMethod2d:
    def __init__(self, call, stepper, title=None, fieldName="pressure"):
        self.module_call = call
        self.module_stepper = stepper
        self.fieldName = fieldName
        if (title == None):
            self.module_title = DummyTitle()
        else:
            self.module_title = title

    def __call__(self, i, j):
        return self.module_call(i, j,self.fieldName)
    
class DummyTitle:
    def __init__(self):
        return
    def __call__(self):
        return "--"
    
class MakeTitle:
    def __init__(self, obj, var, name):
        self.obj = obj
        self.var = var
        self.name = name
    def __call__(self):
        member_value = getattr(self.obj, self.var)
        return f"%s: %03.3e" % (self.name,member_value)