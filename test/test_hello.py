"""
Test hello module for the vspec-template.
"""

import pytest

from vspec_template import hello


def test_hello():
    """
    Test function for the hello() function.
    Asserts that hello() returns the string "Hello, world!".
    """
    assert hello.hello() == "Hello, world!"


def test_add():
    """
    Test function for the add() function.
    Asserts that add(1.0, 2.0) returns 3.0.
    """
    assert hello.add(1.0, 2.0) == pytest.approx(3.0, abs=1e-6)
    with pytest.raises(TypeError):
        hello.add(1.0, "2.0")
    with pytest.raises(TypeError):
        hello.add(1, 2.0)
