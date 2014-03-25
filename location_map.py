#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
from os import system, path
import urllib2, urllib, json
import math
import time

def physical_separation(current_lat, current_long, test_lat, test_long):
    # Inputs are given in degrees, and converted to radians in coordinate transformations
    R = 6378.1*1000 # meters
    current_x = R*math.cos(current_lat*0.01745329252)*math.cos(current_long*0.01745329252)
    current_y = R*math.cos(current_lat*0.01745329252)*math.sin(current_long*0.01745329252)
    current_z = R*math.sin(current_lat*0.01745329252)
    test_x = R*math.cos(test_lat*0.01745329252)*math.cos(test_long*0.01745329252)
    test_y = R*math.cos(test_lat*0.01745329252)*math.sin(test_long*0.01745329252)
    test_z = R*math.sin(test_lat*0.01745329252)
    
    d = math.sqrt( (current_x-test_x)**2 + (current_y-test_y)**2 + (current_z-test_z)**2 )
    # This is a straight-line distance, as if we burrowed through the Earth.
    
    alpha = 2*math.asin(d/(2*R))
    D = R*alpha
    # This the the distance along the surface of the Earth.
    return D

# system("/usr/bin/CoreLocationCLI -once > /tmp/location.txt")
# system("/usr/bin/CoreLocationCLI -once > /tmp/location.txt")
# system("/usr/bin/CoreLocationCLI -once > /tmp/location.txt")

location_file = file("/tmp/location.txt", "r")
line = location_file.readline()
location_file.close()
latitude = line.split("<")[1].split(",")[0].strip("+")
longitude = line.split(">")[0].split(",")[1].strip("+")
sep_meters = 999

modified_time = path.getmtime("/Users/cklein/Desktop/Personal/Misc_Other/coordinates.txt")
hours_ago = (time.time() - modified_time) / 3600.0

if hours_ago < 0.5:
    try:
        old_coord_file = file("/Users/cklein/Desktop/Personal/Misc_Other/coordinates.txt", "r")
        old_latitude_str, old_longitude_str, old_sep_meters = old_coord_file.read().split()
        old_coord_file.close()
        old_latitude = float(old_latitude_str)
        old_longitude = float(old_longitude_str)
        sep_meters = physical_separation(float(latitude), float(longitude), old_latitude, old_longitude)
    except:
        old_coord_file = file("/Users/cklein/Desktop/Personal/Misc_Other/coordinates.txt", "w")
        old_coord_file.write(str(latitude) + "\t" + str(longitude) + "\t" + str(sep_meters) + "\n")
        old_coord_file.close()
        sep_meters = 999
else:
    old_coord_file = file("/Users/cklein/Desktop/Personal/Misc_Other/coordinates.txt", "w")
    old_coord_file.write(str(latitude) + "\t" + str(longitude) + "\t" + str(sep_meters) + "\n")
    old_coord_file.close()
    sep_meters = 999

if sep_meters > 10:
    
    zoom_map_url = """http://maps.googleapis.com/maps/api/staticmap?center=""" + latitude + """,""" + longitude + """&zoom=17&size=512x512&markers=color:blue%7C""" + latitude + """,""" + longitude + """&format=jpg&maptype=roadmap&sensor=false"""
    
    wide_map_url = """http://maps.googleapis.com/maps/api/staticmap?center=""" + latitude + """,""" + longitude + """&zoom=12&size=512x512&markers=color:blue%7C""" + latitude + """,""" + longitude + """&format=jpg&maptype=roadmap&sensor=false"""
        
    system("""/usr/local/bin/wget -q --connect-timeout=1 --output-document=/Users/cklein/Desktop/Personal/Misc_Other/map.jpg \"""" + zoom_map_url + """\"\n""")
    
    system("""/usr/local/bin/wget -q --connect-timeout=1 --output-document=/Users/cklein/Desktop/Personal/Misc_Other/wide_map.jpg \"""" + wide_map_url + """\"\n""")
    
    old_coord_file = file("/Users/cklein/Desktop/Personal/Misc_Other/coordinates.txt", "w")
    old_coord_file.write(str(latitude) + "\t" + str(longitude) + "\t" + str(sep_meters) + "\n")
    old_coord_file.close()
    
    # print str(round(float(latitude), 5)) + "\n" + str(round(float(longitude), 5))
    # print address