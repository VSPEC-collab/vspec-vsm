"""
Miscalaneous helper functions
"""

def round_teff(teff):
    """
    Round the effective temperature to the nearest integer.
    The goal is to reduce the number of unique effective temperatures
    while not affecting the accuracy of the model.

    Parameters
    ----------
    teff : astropy.units.Quantity
        The temperature to round.

    Returns
    -------
    astropy.units.Quantity
        The rounded temperature.

    Notes
    -----
    This function rounds the given effective temperature to the nearest integer value. It is designed to decrease the number of unique effective temperatures while maintaining the accuracy of the model.

    Examples
    --------
    >>> teff = 1234.56 * u.K
    >>> rounded_teff = round_teff(teff)
    >>> print(rounded_teff)
    1235 K

    >>> teff = 2000.4 * u.K
    >>> rounded_teff = round_teff(teff)
    >>> print(rounded_teff)
    2000 K

    """

    val = teff.value
    unit = teff.unit
    return int(round(val)) * unit