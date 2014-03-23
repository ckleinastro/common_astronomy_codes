#!/Users/cklein/anaconda/bin/python
import matplotlib 
matplotlib.use('Agg')

from os import system
import ephem
import operator

from scipy import *
import pylab
from pylab import *
# from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

from math import radians
from numpy import pi, arctan2, sin, cos, arcsin

# system("wget http://sv.berkeley.edu:/view/images/newview.jpg .")

system("/usr/bin/CoreLocationCLI -once > /tmp/location.txt")
system("/usr/bin/CoreLocationCLI -once > /tmp/location.txt")
system("/usr/bin/CoreLocationCLI -once > /tmp/location.txt")
location_file = file("/tmp/location.txt", "r")
line = location_file.readline()
location_file.close()
latitude = line.split("<")[1].split(",")[0].strip("+")
longitude = line.split(">")[0].split(",")[1]

observatory = ephem.Observer()
observatory.long = longitude
observatory.lat = latitude
# Berkeley Coordinates
# observatory.long = "-122.25803505"
# observatory.lat = "37.86934626"
observatory.elevation = 200

right_now = ephem.now()

observatory.date = right_now

sun = ephem.Sun()
sun.compute(observatory)
moon = ephem.Moon()
moon.compute(observatory)

# days since last full moon
last_full_moon_time = ephem.previous_full_moon(right_now) - right_now

# days until next full moon
next_full_moon_time = ephem.next_full_moon(right_now) - right_now

# observatory.date = right_now

events_list = []

sunset_time = observatory.next_setting(sun)
observatory.date = right_now
events_list.append(["Sunset: \t", ephem.localtime(sunset_time).strftime("%I:%M%p"), sunset_time - right_now])

sunnoon_time = observatory.next_transit(sun)
observatory.date = right_now
events_list.append(["Local Noon:", ephem.localtime(sunnoon_time).strftime("%I:%M%p"), sunnoon_time - right_now])

sunrise_time = observatory.next_rising(sun)
observatory.date = right_now
events_list.append(["Sunrise:\t", ephem.localtime(sunrise_time).strftime("%I:%M%p"), sunrise_time - right_now])

moonset_time = observatory.next_setting(moon)
observatory.date = right_now
events_list.append(["Moonset:", ephem.localtime(moonset_time).strftime("%I:%M%p"), moonset_time - right_now])

moonnoon_time = observatory.next_transit(moon)
observatory.date = right_now
events_list.append(["Moon Noon:", ephem.localtime(moonnoon_time).strftime("%I:%M%p"), moonnoon_time - right_now])

moonrise_time = observatory.next_rising(moon)
observatory.date = right_now
events_list.append(["Moonrise:", ephem.localtime(moonrise_time).strftime("%I:%M%p"), moonrise_time - right_now])

events_list_sorted = sorted(events_list, key=operator.itemgetter(2))

system("/Users/cklein/anaconda/bin/python /Users/cklein/Desktop/Personal/Misc_Other/location_map.py")

# system("/Users/cklein/anaconda/bin/python /Users/cklein/Desktop/Personal/Misc_Other/constellation_drawer/draw_constellation.py")

# events_file = file("/Users/cklein/Desktop/Personal/Misc_Other/celestial_timings.txt", "w")

for entry in events_list_sorted:
    print entry[0] + "\t" + entry[1].lstrip("0")
    # events_file.write(entry[0] + "\t" + entry[1].lstrip("0") + "\n")
print "Prev Full Moon: " + str(round(last_full_moon_time, 2)) + "d"
# events_file.write("Prev Full Moon: " + str(round(last_full_moon_time, 2)) + "d" + "\n")
print "Next Full Moon: +" + str(round(next_full_moon_time, 2))+ "d"
# events_file.write("Next Full Moon: +" + str(round(next_full_moon_time, 2))+ "d" + "\n")
# events_file.close()

# Before running Messages status updater, make sure that Messages is actually running.
system("ps aux | grep Messages > /tmp/Messages_process.txt")
process_file = file("/tmp/Messages_process.txt", "r")
for line in process_file:
    if line.split()[10] == "/Applications/Messages.app/Contents/MacOS/Messages":
        system("/Users/cklein/anaconda/bin/python /Users/cklein/Desktop/Personal/Misc_Other/Messages_status_update.py")
process_file.close()


def gal2eq(l, b):
    # convert first from degrees to radians
    l = radians(l)
    b = radians(b)
    
    # Here we convert from Galactic to Equatorial (all conversions done in radians)
    # North galactic pole (J2000) -- according to Wikipedia
    pole_ra = radians(192.859508)
    pole_dec = radians(27.128336)
    posangle = radians(122.932-90.0)
    ra = arctan2( (cos(b)*cos(l-posangle)), (sin(b)*cos(pole_dec) - cos(b)*sin(pole_dec)*sin(l-posangle)) ) + pole_ra
    dec = arcsin( cos(b)*cos(pole_dec)*sin(l-posangle) + sin(b)*sin(pole_dec) )

    # Convert from radians to degrees for output in degrees
    ra = ra * 180./pi # degrees
    dec = dec * 180./pi # degrees

    return (ra, dec)



limiting_mag = 5.0


observatory.date = ephem.now() + 3./60./24 # add three minutes to account for update lag


new_point = ephem.Equatorial(str(12), str(90), epoch="2000")
ephem_line = ("SouthPole" + ",f|M|F7," + str(new_point.ra) + "," + 
    str(new_point.dec) + ",0,2000")
star_calced = ephem.readdb(ephem_line)
star_calced.compute(observatory)
north_pole_alt = star_calced.alt
north_pole_az = star_calced.az
north_pole_zenith_radii = 90 - north_pole_alt*180/pi
north_pole_az = north_pole_az + pi/2
north_pole_plot_x = north_pole_zenith_radii * cos(north_pole_az)
north_pole_plot_y = north_pole_zenith_radii * sin(north_pole_az)

new_point = ephem.Equatorial(str(12), str(-90), epoch="2000")
ephem_line = ("SouthPole" + ",f|M|F7," + str(new_point.ra) + "," + 
    str(new_point.dec) + ",0,2000")
star_calced = ephem.readdb(ephem_line)
star_calced.compute(observatory)
south_pole_alt = star_calced.alt
south_pole_az = star_calced.az
south_pole_zenith_radii = 90 - south_pole_alt*180/pi
south_pole_az = south_pole_az + pi/2
south_pole_plot_x = south_pole_zenith_radii * cos(south_pole_az)
south_pole_plot_y = south_pole_zenith_radii * sin(south_pole_az)


star_data = loadtxt("/Users/cklein/Desktop/Personal/Misc_Other/star_catalog.txt")

star_alts = []
star_azs = []
star_mags = []
for star in star_data:
    if star[3] < limiting_mag:
        new_point = ephem.Equatorial(str(star[1]/15.), str(star[2]), epoch="2000")
        ephem_line = ("SAO"+str(star[0]) + ",f|M|F7," + str(new_point.ra) + "," + 
            str(new_point.dec) + ",0,2000")
        star_calced = ephem.readdb(ephem_line)
        star_calced.compute(observatory)
        if star_calced.alt > 0:
            star_alts.append(star_calced.alt)
            star_azs.append(star_calced.az)
            star_mags.append(star[3])

star_azs = array(star_azs) # similar to angle around the circle
star_alts = array(star_alts) # similar to radius from center (zenith)
star_zenith_radii = 90 - star_alts*180/pi
star_azs = star_azs + pi/2
star_plot_x = star_zenith_radii * cos(star_azs)
star_plot_y = star_zenith_radii * sin(star_azs)




equator_alts = []
equator_azs = []

equator_ras = linspace(0,360,1200)
equator_decs = zeros(1200)

for n in range(len(equator_ras)):
    new_point = ephem.Equatorial(str(equator_ras[n]/15.), str(equator_decs[n]), epoch="2000")
    ephem_line = ("EquatorSeg"+",f|M|F7," + str(new_point.ra) + "," + 
        str(new_point.dec) + ",0,2000")
    equator_calced = ephem.readdb(ephem_line)
    equator_calced.compute(observatory)
    if equator_calced.alt > 0:
        equator_alts.append(equator_calced.alt)
        equator_azs.append(equator_calced.az)

equator_azs = array(equator_azs) # similar to angle around the circle
equator_alts = array(equator_alts) # similar to radius from center (zenith)
equator_zenith_radii = 90 - equator_alts*180/pi
equator_azs = equator_azs + pi/2
equator_plot_x = equator_zenith_radii * cos(equator_azs)
equator_plot_y = equator_zenith_radii * sin(equator_azs)

equator_plot_data = zip(equator_plot_x, equator_plot_y)
equator_plot_data.sort()
equator_plot_data = array(equator_plot_data)

# Create the line of points that represent the Galactic Equator
def make_gal_equator_line(b=0, sigma=60):
    try:
        gal_l = linspace(-180,180,400)
        gal_equator_ras = []
        gal_equator_decs = []
        for l in gal_l:
            gal_eq_ra, gal_eq_dec = gal2eq(l, ((b*20/(sigma*(2*pi)**0.5))*exp(-0.5*(l/sigma)**2)))
            gal_equator_ras.append(gal_eq_ra)
            gal_equator_decs.append(gal_eq_dec)
        gal_equator_alts = []
        gal_equator_azs = []
        for n in range(len(gal_equator_ras)):
            new_point = ephem.Equatorial(str(gal_equator_ras[n]/15.), str(gal_equator_decs[n]), epoch="2000")
            ephem_line = ("EquatorSeg"+",f|M|F7," + str(new_point.ra) + "," + 
                str(new_point.dec) + ",0,2000")
            gal_equator_calced = ephem.readdb(ephem_line)
            gal_equator_calced.compute(observatory)
            if gal_equator_calced.alt > 0:
                gal_equator_alts.append(gal_equator_calced.alt)
                gal_equator_azs.append(gal_equator_calced.az)
        gal_equator_azs = array(gal_equator_azs) # similar to angle around the circle
        gal_equator_alts = array(gal_equator_alts) # similar to radius from center (zenith)
        gal_equator_zenith_radii = 90 - gal_equator_alts*180/pi
        gal_equator_azs = gal_equator_azs + pi/2
        gal_equator_plot_x = gal_equator_zenith_radii * cos(gal_equator_azs)
        gal_equator_plot_y = gal_equator_zenith_radii * sin(gal_equator_azs)
        gal_equator_plot_data = zip(gal_equator_plot_x, gal_equator_plot_y)
        # gal_equator_plot_data.sort()
        gal_equator_plot_data = array(gal_equator_plot_data)
        previous_point = array([])
        sep_list_x = []
        sep_list_y = []
        for point in gal_equator_plot_data:
            if previous_point.any():
                sep_list_x.append(point[0] - previous_point[0])
                sep_list_y.append(point[1] - previous_point[1])
            previous_point = point
        sep_array_x = array(sep_list_x)
        sep_array_y = array(sep_list_y)
        sep_array = (sep_array_x**2 + sep_array_y**2)**0.5
        if max(abs(sep_array)) > 2:
            while max(abs(sep_array)) > 2:
                jump_index = where(abs(sep_array) == max(abs(sep_array)))[0][0]+1
                gal_equator_plot_data_reordered = append(gal_equator_plot_data[jump_index:], gal_equator_plot_data[:jump_index], axis=0)
                gal_equator_plot_data = gal_equator_plot_data_reordered
                previous_point = array([])
                sep_list_x = []
                sep_list_y = []
                for point in gal_equator_plot_data:
                    if previous_point.any():
                        sep_list_x.append(point[0] - previous_point[0])
                        sep_list_y.append(point[1] - previous_point[1])
                    previous_point = point
                sep_array_x = array(sep_list_x)
                sep_array_y = array(sep_list_y)
                sep_array = (sep_array_x**2 + sep_array_y**2)**0.5
            
        return gal_equator_plot_data
    except:
        return array([[]])
        
gal_equator_plot_data = make_gal_equator_line(b=0)[1:-1]
gal_equator_plot_data_plus5 = make_gal_equator_line(b=50)[1:-1]
gal_equator_plot_data_minus5 = make_gal_equator_line(b=-50)[1:-1]

gal_center_visible = False
gal_center_ras = []
gal_center_decs = []

gal_eq_radius = 3
l_circ = gal_eq_radius*sin(linspace(0,2*pi,100))
b_circ = gal_eq_radius*cos(linspace(0,2*pi,100))
for n in range(len(l_circ)):
    gal_center_ra, gal_center_dec = gal2eq(l_circ[n], b_circ[n])
    gal_center_ras.append(gal_center_ra)
    gal_center_decs.append(gal_center_dec)

gal_center_alts = []
gal_center_azs = []
for n in range(len(gal_center_ras)):
    new_point = ephem.Equatorial(str(gal_center_ras[n]/15.), str(gal_center_decs[n]), epoch="2000")
    ephem_line = ("centerSeg"+",f|M|F7," + str(new_point.ra) + "," + 
        str(new_point.dec) + ",0,2000")
    gal_center_calced = ephem.readdb(ephem_line)
    gal_center_calced.compute(observatory)
    if gal_center_calced.alt*180/pi > gal_eq_radius:
        gal_center_alts.append(gal_center_calced.alt)
        gal_center_azs.append(gal_center_calced.az)
if len(gal_center_alts) == 100:
    gal_center_visible = True
gal_center_azs = array(gal_center_azs) # similar to angle around the circle
gal_center_alts = array(gal_center_alts) # similar to radius from center (zenith)
gal_center_zenith_radii = 90 - gal_center_alts*180/pi
gal_center_azs = gal_center_azs + pi/2
gal_center_plot_x = gal_center_zenith_radii * cos(gal_center_azs)
gal_center_plot_y = gal_center_zenith_radii * sin(gal_center_azs)
gal_center_plot_data = zip(gal_center_plot_x, gal_center_plot_y)
gal_center_plot_data = array(gal_center_plot_data)


remove_stars = [
    "Alcyone", 
    "Merope",
    "Maia",
    "Taygeta",
    "Electra",
    "Bellatrix",
    "Alnilam",
    "Alnitak",
    "Wezen",
    "Arkab Posterior",
    ]

starnames = []
named_star_alts = []
named_star_azs = []
import ephem.stars 
for star in ephem.stars.db.split("\n"):
    starname = star.split(",")[0]
    if len(starname)>0 and starname not in remove_stars:
        star_calced = ephem.star(starname)
        star_calced.compute(observatory)
        if star_calced.alt > 0:
            starnames.append(starname)
            named_star_alts.append(star_calced.alt)
            named_star_azs.append(star_calced.az)
named_star_azs = array(named_star_azs) # similar to angle around the circle
named_star_alts = array(named_star_alts) # similar to radius from center (zenith)
named_star_zenith_radii = 90 - named_star_alts*180/pi
named_star_azs = named_star_azs + pi/2
named_star_plot_x = named_star_zenith_radii * cos(named_star_azs)
named_star_plot_y = named_star_zenith_radii * sin(named_star_azs)


planets = ["Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Moon", "Sun"]
planet_names = []
planet_alts = []
planet_azs = []
for planet in planets:
    exec "planet_calced = ephem." + planet + "()"
    planet_calced.compute(observatory)
    if planet_calced.alt > 0:
        planet_names.append(planet)
        planet_alts.append(planet_calced.alt)
        planet_azs.append(planet_calced.az)
planet_azs = array(planet_azs) # similar to angle around the circle
planet_alts = array(planet_alts) # similar to radius from center (zenith)
planet_zenith_radii = 90 - planet_alts*180/pi
planet_azs = planet_azs + pi/2
planet_plot_x = planet_zenith_radii * cos(planet_azs)
planet_plot_y = planet_zenith_radii * sin(planet_azs)



zenith_ra, zenith_dec = observatory.radec_of(0,pi/2)

close("all")

fig = plt.figure(figsize=(12, 12))





ax = subplot(111, aspect="equal")

# Plot the lines of constant latitude
circle_azs = linspace(0,2*pi,1000)
lat_radius_list = [90, 60, 30]
for lat_radius in lat_radius_list:
    horizon_line_x = lat_radius*cos(circle_azs)
    horizon_line_y = lat_radius*sin(circle_azs)
    ax.plot(horizon_line_x, horizon_line_y, linestyle="--", color="gray", alpha=1, linewidth=1)
#     if lat_radius < 90:
#         ax.text(lat_radius/sqrt(2), lat_radius/sqrt(2), str(90-lat_radius), alpha=0.5,
#             horizontalalignment="left", verticalalignment="bottom", color="white")
#         ax.text(lat_radius/sqrt(2), -lat_radius/sqrt(2), str(90-lat_radius), alpha=0.5, 
#             horizontalalignment="left", verticalalignment="top", color="white")
#         ax.text(-lat_radius/sqrt(2), lat_radius/sqrt(2), str(90-lat_radius), alpha=0.5,
#             horizontalalignment="right", verticalalignment="bottom", color="white")
#         ax.text(-lat_radius/sqrt(2), -lat_radius/sqrt(2), str(90-lat_radius), alpha=0.5,
#             horizontalalignment="right", verticalalignment="top", color="white")


ax.plot(equator_plot_data[:,0], equator_plot_data[:,1], color="green", alpha=0.7)

try:
    ax.plot(gal_equator_plot_data[:,0], gal_equator_plot_data[:,1], color="HotPink", alpha=0.7)
except:
    line_fail=True
try:
    ax.plot(gal_equator_plot_data_plus5[:,0], gal_equator_plot_data_plus5[:,1], color="HotPink", alpha=0.7, linestyle="--")
except:
    line_fail=True
try:
    ax.plot(gal_equator_plot_data_minus5[:,0], gal_equator_plot_data_minus5[:,1], color="HotPink", alpha=0.7, linestyle="--")
except:
    line_fail=True


if gal_center_visible:
    ax.plot(gal_center_plot_data[:,0], gal_center_plot_data[:,1], color="red", alpha=0.7)

sizes = ((7-array(star_mags))**3)/5
ax.scatter(star_plot_x, star_plot_y, marker="o", color="white", s=sizes)
# ax.scatter(star_plot_x, star_plot_y, marker="o", color="black")

ax.scatter(named_star_plot_x, named_star_plot_y, color="white", marker="o", alpha=1, s=56)
ax.scatter(named_star_plot_x, named_star_plot_y, color="lime", marker="o", alpha=1, s=25)
for n in range(len(starnames)):
    if named_star_plot_x[n] < 75 or abs(named_star_plot_y[n]) > 50:
        ax.text(named_star_plot_x[n]+1.75, named_star_plot_y[n]-1.75, starnames[n], color="white")

ax.scatter(planet_plot_x, planet_plot_y, color="cyan", marker="o", alpha=1, s=50)
for n in range(len(planet_names)):
    if planet_names[n] == "Sun":
        ax.scatter(planet_plot_x[n], planet_plot_y[n], color="orange", marker="o", alpha=1, s=256)
        ax.text(planet_plot_x[n]+2.75, planet_plot_y[n]-2.75, planet_names[n], color="orange")
    else:
        ax.text(planet_plot_x[n]+1.75, planet_plot_y[n]-1.75, planet_names[n], color="cyan")


ax.set_xlim(-93, 93)
ax.set_ylim(-93, 93)

ax.text(-93, 0, "E", horizontalalignment="center", verticalalignment="center", color="white")
ax.text(93, 0, "W", horizontalalignment="center", verticalalignment="center", color="white")
ax.text(0, -93, "S", horizontalalignment="center", verticalalignment="center", color="white")
ax.text(0, 93, "N", horizontalalignment="center", verticalalignment="center", color="white")

ax.text(-90, 90, "Zenith Position", horizontalalignment="left", verticalalignment="center", color="white")
ax.text(-90, 86, r'$\alpha$=' + str(zenith_ra), horizontalalignment="left", verticalalignment="center", color="white")
ax.text(-90, 82, r'$\delta$=' + str(zenith_dec), horizontalalignment="left", verticalalignment="center", color="white")

ax.text(65, 90, "Zenith Position", horizontalalignment="left", verticalalignment="center", color="white")
ax.text(65, 86, r'$\alpha$=' + str(180/pi*zenith_ra)[:10], horizontalalignment="left", verticalalignment="center", color="white")
ax.text(65, 82, r'$\delta$=' + str(180/pi*zenith_dec)[:10], horizontalalignment="left", verticalalignment="center", color="white")

ax.text(-90, -85, "Compute Time", horizontalalignment="left", verticalalignment="center", color="white")
ax.text(-90, -90, str(observatory.date) + " UT", horizontalalignment="left", verticalalignment="center", color="white")
local_date = ephem.localtime(observatory.date)
local_date_str = "%d/%d/%d %02d:%02d:%02d local" % ( local_date.year, local_date.month, local_date.day, local_date.hour, local_date.minute, int(round(local_date.second)) )
ax.text(-90, -95, local_date_str, horizontalalignment="left", verticalalignment="center", color="white")

ax.text(50, -85, "Compute Long & Lat", horizontalalignment="left", verticalalignment="center", color="white")
ax.text(50, -90, str(observatory.long) + " " + str(observatory.lat), horizontalalignment="left", verticalalignment="center", color="white")
ax.text(50, -95, str(observatory.long*180/pi)[:10] + " " + str(observatory.lat*180/pi)[:10], horizontalalignment="left", verticalalignment="center", color="white")

ax.scatter(0,0, color="red", marker="+")
ax.scatter(north_pole_plot_x, north_pole_plot_y, color="magenta", marker='^')
ax.text(north_pole_plot_x, north_pole_plot_y+2, "NP", horizontalalignment="center", verticalalignment="center", color="magenta")
ax.scatter(south_pole_plot_x, south_pole_plot_y, color="magenta", marker='v')
ax.text(south_pole_plot_x, south_pole_plot_y-2, "SP", horizontalalignment="center", verticalalignment="center", color="magenta")


ax.set_xticks([])
ax.set_yticks([])
ax.tick_params(axis='x', labelbottom=False, labeltop=False)
ax.tick_params(axis='y', labelleft=False, labelright=False)

ax.patch.set_facecolor('None')

axis('off')
# fig.set_facecolor("#000000")
# canvas = FigureCanvas(fig)
# canvas.figure.set_facecolor("#000000")
# canvas.print_figure("/Users/cklein/Desktop/Personal/Misc_Other/starmap.pdf", transparent=True)


savefig("/Users/cklein/Desktop/Personal/Misc_Other/starmap.pdf", transparent=True)
close("all")
