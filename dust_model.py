from scipy import *
from scipy.interpolate import interp1d
from scipy import random, stats, integrate, inf
from pylab import *
import sys
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas


# For Bulge stars, we give them an exponential distribution about the GC
def make_bulge(n_stars_bulge, exponential_bulge_parameter=1500):
    pos_gc_bulge_phi = uniform(low=0, high=2*pi, size=n_stars_bulge)
    pos_gc_bulge_costheta = uniform(low=-1, high=1, size=n_stars_bulge)
    pos_gc_bulge_u = uniform(low=0, high=1, size=n_stars_bulge)
    pos_gc_bulge_theta = arccos( pos_gc_bulge_costheta )
    pos_gc_bulge_r = exponential(scale=1.33*exponential_bulge_parameter, size=n_stars_bulge) *  pos_gc_bulge_u**(1/3.)
    pos_gc_bulge_x = pos_gc_bulge_r * sin( pos_gc_bulge_theta) * cos( pos_gc_bulge_phi )
    pos_gc_bulge_y = pos_gc_bulge_r * sin( pos_gc_bulge_theta) * sin( pos_gc_bulge_phi )
    pos_gc_bulge_z = pos_gc_bulge_r * cos( pos_gc_bulge_theta )
    return pos_gc_bulge_x, pos_gc_bulge_y, pos_gc_bulge_z



# For Disk stars we define the position originally in cylindrical coordinates, 
# with exponential radial, r, and height, z, positions.
def make_disk(n_stars_disk, exponential_disk_r_parameter=3000, exponential_disk_z_parameter=250):
    pos_gc_disk_phi = uniform(low=0, high=2*pi, size=n_stars_disk)
    pos_gc_disk_r = exponential(scale=exponential_disk_r_parameter, size=n_stars_disk) * where(uniform(size=n_stars_disk) < 0.5, -1, 1)
    pos_gc_disk_z = exponential(scale=exponential_disk_z_parameter, size=n_stars_disk) * where(uniform(size=n_stars_disk) < 0.5, -1, 1)
    pos_gc_disk_x = pos_gc_disk_r * cos(pos_gc_disk_phi)
    pos_gc_disk_y = pos_gc_disk_r * sin(pos_gc_disk_phi)
    return pos_gc_disk_x, pos_gc_disk_y, pos_gc_disk_z



exponential_distribution = lambda l, x: exp(-x/l)/l

gaussian_distribution = lambda s, x: (1/(s*sqrt(2*pi))) * exp(-1 * (x**2)/(2*s**2))



galaxy_radial_profile = lambda a, s, b, l, c, x: a*( (b/(s*sqrt(2*pi))) * exp(-1 * (x**2)/(2*s**2)) + c*exp(-x/l)/l )


x = linspace(-15000, 15000, 1000)
y = linspace(-15000, 15000, 1000)
z = linspace(-15000, 15000, 1000)


# Dust density is gaussian bulge in spherical radius (sqrt(x**2+y**2+z**2)), exponential disk in cylindrical radius (sqrt(x**2+y**2)), and gaussian in height (z).
#   a   normalization constant for the full distribution
#   s   width (sigma) of the spherical gaussian bulge
#   b   height (or strength) of the spherical gaussian bulge
#   l   scale length of the exponential disk
#   c   height (or strength) of the exponential disk
#   d   height (or strength) of the gaussian height (z) profile
#   e   width (sigma) of the gaussian height (z) profile
dust_density = lambda a, s, b, l, c, d, e, z, y, x: a*( (b/(s*sqrt(2*pi))) * exp(-1 * (sqrt(x**2+y**2+z**2)**2)/(2*s**2)) + c*exp(-sqrt(x**2+y**2)/l)/l ) * (d/(e*sqrt(2*pi))) * exp(-1 * (z**2)/(2*e**2))

dd = dust_density(672.00089811262512, 1000., 5., 4000., 7., 300., 300., z, y, x)

dust = lambda z, y, x: (0.004320427207007109)*( (5./(1000.*sqrt(2.*pi))) * exp(-1. * (sqrt(x**2.+y**2.+z**2.)**2.)/(2.*1000.**2.)) + 7.*exp(-sqrt(x**2.+y**2)/4000.)/4000.) * (300./(300.*sqrt(2*pi))) * exp(-1. * (z**2.)/(2.*300.**2.))


sun_gc_x = 8330.
sun_gc_y = 0.
sun_gc_z = 25.

test_x = 1.
test_y = 1.
test_z = 1.


"""
New notes for improving the dust integration.
80 pc in exponential disk height
Bulge seems too large (nicholiave and weinberg used 2mass data) Mike Rich from UCLA.

Also do disk scale height and length recovery (using only disk populations) with
new simulations for all 3 survey durations. Use dust, but also use the disk pizza
slicing. ("plot_WITS.py" after running the simulation)

One galaxy where the 3 PV classes have the same disk scales, and one where they
are different. Convince ourselves that we can distinguish between those two.

"""

def integrate_dust(test_x, test_y, test_z, sun_gc_x=8330., sun_gc_y=0., sun_gc_z=25.):
    # Dust density is gaussian bulge in spherical radius (sqrt(x**2+y**2+z**2)), exponential disk in cylindrical radius (sqrt(x**2+y**2)), and gaussian in height (z).
    a=0.00011672834801957508 #  normalization constant for the full distribution. Fixed so that line of sight extinction from the GC=1.5.
    s=500.                  #  width (sigma, pc) of the spherical gaussian bulge
    b=5.                    #  height (or strength) of the spherical gaussian bulge
    l=4000.                 #  scale length (pc) of the exponential disk
    c=7.                    #  height (or strength) of the exponential disk
    d=300.                  #  height (or strength) of the gaussian height (z) profile
    e=80.                   #  width (sigma, pc) of the gaussian height (z) profile
    
    dust = lambda z, y, x: a*( (b/(s*sqrt(2*pi))) * exp(-1 * (sqrt(x**2+y**2+z**2)**2)/(2*s**2)) + c*exp(-sqrt(x**2+y**2)/l)/l ) * (d/(e*sqrt(2*pi))) * exp(-1 * (z**2)/(2*e**2))
    
    # dust = lambda z, y, x: (0.0011589713316126277)*( (5./(1000.*sqrt(2.*pi))) * exp(-1. * (sqrt(x**2.+y**2.+z**2.)**2.)/(2.*1000.**2.)) + 7.*exp(-sqrt(x**2.+y**2)/4000.)/4000.) * (300./(80.*sqrt(2*pi))) * exp(-1. * (z**2.)/(2.*80.**2.))
    
    x_array = linspace(test_x, sun_gc_x, 1000)
    y_array = linspace(test_y, sun_gc_y, 1000)
    z_array = linspace(test_z, sun_gc_z, 1000)
    a = dust(z_array, y_array, x_array) * sqrt( (test_x - sun_gc_x)**2 + (test_y - sun_gc_y)**2 + (test_z - sun_gc_z)**2 )
    a_err = max(0.1*a.sum(), 0.01)
    return a.sum(), a_err


num_pix = 101
map_extent=20000

x_distances = linspace(-map_extent, map_extent, num_pix)
y_distances = linspace(-map_extent, map_extent, num_pix)



rcParams['axes.unicode_minus'] = False
matplotlib.rc('font', family="serif")
matplotlib.rc('font', serif="Times New Roman")
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rc('axes', labelsize=16)


z_val = 0.0

for z_val in linspace(-500, 500, 41):


    extinction_list = []
    for y_d in y_distances:
        for x_d in x_distances:
            los_extinction, los_extinction_errs = integrate_dust(x_d, y_d, z_val)
            extinction_list.append(los_extinction)

    extinction_array = array(extinction_list).reshape((num_pix,num_pix))

    """ Here we make the color mapping to color the distances """

    # This data was extracted from Carl's plot of his colormap
    raw_red_color = loadtxt("/Users/cklein/Desktop/Research_Support/python_scripts/color_mapping_code/red_mapping.txt")
    raw_green_color = loadtxt("/Users/cklein/Desktop/Research_Support/python_scripts/color_mapping_code/green_mapping.txt")
    raw_blue_color = loadtxt("/Users/cklein/Desktop/Research_Support/python_scripts/color_mapping_code/blue_mapping.txt")

    red_function = interp1d(raw_red_color[:,0], raw_red_color[:,1], kind="linear")
    green_function = interp1d(raw_green_color[:,0], raw_green_color[:,1], kind="linear")
    blue_function = interp1d(raw_blue_color[:,0], raw_blue_color[:,1], kind="linear")

    # now create the 2D colormap, multiplying by blackness to dimm vertically
    color_display = zeros((256, 256, 3))
    color_display[:,:,0] = (red_function(arange(256*256).reshape(256,256)/256).transpose()  / 255.)
    color_display[:,:,1] = (green_function(arange(256*256).reshape(256,256)/256).transpose()  / 255.)
    color_display[:,:,2] = (blue_function(arange(256*256).reshape(256,256)/256).transpose()  / 255.)

    def make_color(val):
        val *= 255
        color = [red_function(val)/255., green_function(val)/255., blue_function(val)/255., 1.0]
        return color

    extinction_color_vals = (array(extinction_list) - 0.0) / (3.1- 0.0)

    extinction_colors_list = []
    for val in extinction_color_vals:
        extinction_colors_list.append(make_color(val))

    extinction_colors_array = array(extinction_colors_list).reshape((num_pix,num_pix, 4))



    pc_per_pixel = 2.*map_extent / float(num_pix-1.)
    gc_x_pix = (num_pix-1.)/2.
    gc_y_pix = (num_pix-1.)/2.
    sun_x_pix = gc_x_pix + 8330./pc_per_pixel
    sun_y_pix = gc_y_pix

    fig = plt.figure(figsize=(7, 7))
    ax1 = fig.add_subplot(2,1,1)

    ax1.imshow(extinction_colors_array, origin="lower", interpolation="lanczos")

    ax1.set_xlabel("X [pc]")
    ax1.set_ylabel("Y [pc]")
    ax1.set_title("WITS Line of Sight Extinction at midplane height Z=" + str(int(z_val)))
    major_tick_spacing_pixels = 25
    majorLocator_x = MultipleLocator(major_tick_spacing_pixels)
    minorLocator_x = MultipleLocator(5)
    ax1.xaxis.set_major_locator(majorLocator_x)
    ax1.xaxis.set_minor_locator(minorLocator_x)

    majorLocator_y1 = MultipleLocator(major_tick_spacing_pixels)
    minorLocator_y1 = MultipleLocator(10)
    ax1.yaxis.set_major_locator(majorLocator_y1)
    ax1.yaxis.set_minor_locator(minorLocator_y1)

    float_tick_labels = arange(0, num_pix, major_tick_spacing_pixels)*pc_per_pixel - map_extent
    str_tick_labels = []
    for entry in float_tick_labels:
        str_tick_labels.append(str(int(entry)))
    ax1.xaxis.set_ticklabels(str_tick_labels)
    ax1.yaxis.set_ticklabels(str_tick_labels)

    ax1.scatter(sun_x_pix, sun_y_pix, color="white", marker=(5,1), s=60)
    ax1.scatter(sun_x_pix, sun_y_pix, color="blue", marker=(5,1), s=12)
    ax1.scatter(gc_x_pix, gc_y_pix, color="white", marker=(4,1), s=60)
    ax1.scatter(gc_x_pix, gc_y_pix, color="red", marker=(4,1), s=12)



    ax1.set_xlim(-0.5, num_pix-0.5)
    ax1.set_ylim(-0.5, num_pix-0.5)


    (array(extinction_list) - 0.0) / (3.1- 0.0)

    ax2 = fig.add_subplot(2,1,2)
    ax2.imshow(color_display, origin="lower", interpolation="lanczos", 
        extent=[0.0, 3.1, 0, 
        0.05*(3.1- 0.0)], alpha=1)
    ax2.set_xlabel(r"$A_{\bf{W1}}$")
    ax2.xaxis.set_label_position("bottom")
    ax2.yaxis.set_ticks([])
    # [left, bottom, width, height]
    ax1.set_position([0.155, 0.175, 0.75, 0.8])
    ax2.set_position([0.155, 0.045, 0.75, 0.10])

    canvas = FigureCanvas(fig)
    canvas.print_figure("extinction_plots/los_extinction_z" + str(z_val) + "_new.jpg", dpi=144)
    close("all")
