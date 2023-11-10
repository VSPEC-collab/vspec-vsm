"""
Hello module for the vspec-template
-----------------------------------

This is the hello module for the vspec-template.
It does not do much.

For info on writing numpydoc style docstrings, see:
https://numpydoc.readthedocs.io/en/latest/format.html
"""


def _hello():
    """
    Private function to say "Hello, world!"

    Private functions like this do not show up in the documentation,
    but they can still be documented for internal use.

    Returns
    -------
    str
        A string to say "Hello, world!"

    Notes
    -----
    This is an example of 
    numpydoc style docstrings.
    """
    return "Hello, world!"


def hello():
    """
    Hello world! function.


    Returns
    -------
    str
        A string that says "Hello, world!"
    """
    return _hello()


def add(a: float, b: float) -> float:
    """
    Adds two float numbers together and returns the sum.

    Parameters
    ----------
    a : float
        The first number to be added.
    b : float
        The second number to be added.

    Returns
    -------
    float
        The sum of the two input numbers.

    Raises
    ------
    TypeError
        If either `a` or `b` is not of type float.
    """
    if not isinstance(a, float):
        msg = "a must be a float"
        raise TypeError(msg)
    if not isinstance(b, float):
        msg = "b must be a float"
        raise TypeError(msg)
    return a + b
