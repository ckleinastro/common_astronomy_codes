import ephem

# Calculate the heliocentric julian date.
def heliocentric_julian_date(observation_time, observation_ra_radians, 
    observation_dec_radians):
    # Compute the observation time in Heliocentric Julian Date. First convert to 
    # julian date (jd) by adding the offset to ephem.date format.
    observation_time_jd = float(observation_time) + 2415020.0
    # Calculate the Sun's position in the sky.
    sun = ephem.Sun()
    sun.compute(observation_time)
    sun_ra_radians = float(sun.a_ra)
    sun_dec_radians = float(sun.a_dec)
    # Calculate the Earth-Sun light travel time in days.
    earth_sun_light_travel_time = sun.earth_distance * 0.00577551833
    # Calculate the observation time in heliocentric julian date.
    observation_time_hjd = (observation_time_jd - earth_sun_light_travel_time * 
        (sin(observation_dec_radians) * sin(sun_ra_radians) + 
        cos(observation_dec_radians) * cos(sun_ra_radians) * 
        cos(observation_dec_radians - sun_dec_radians)))
    return observation_time_hjd