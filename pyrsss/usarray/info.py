import pandas as PD


INFO_LINK = 'http://ds.iris.edu/files/earthscope/usarray/_US-MT-StationList.txt'
"""Full URL link to USArray MT site information text file."""


HEADER = 'VNET	NET	STA	SITE DESCRIPTION	LAT	LON	ELEV	START	END	STATUS	INSTALL	CERT'
"""Header line of data file."""


def info_map(info_link=INFO_LINK):
    """
    Return a :class:`DataFrame` containing the information provided at
    *info_link*, a link to a tab delineated text file containing
    information for each USArray MT site.
    """
    df = PD.read_table(info_link,
                       sep='\t',
                       skiprows=1,
                       names=['vnet',
                              'net',
                              'sta',
                              'location',
                              'lat',
                              'lon',
                              'elev',
                              'start',
                              'end',
                              'status',
                              'install',
                              'cert'],
                       parse_dates=[7, 8],
                       index_col=2)
    return df
