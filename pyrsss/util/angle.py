def convert_lon(lon):
    """
    Return 0 <= *lon* < 360 converted to -180 <= *lon < 180.
    """
    return lon - 360 if lon > 180 else lon
