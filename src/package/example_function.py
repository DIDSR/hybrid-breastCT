"""
A template script.
"""

import numpy as np
from collections import OrderedDict


def foo(message: str, text: str = "default_value") -> tuple[str, int]:
    """
    Sample docstring.

    Args:
        message (str): text to be processed.
        text (str, optional): _description_. Defaults to "default_value".

    Returns:
        tuple[str, int]: _description_
    """
    output_text = message.upper() + text
    output_text_length = len(output_text)

    return output_text, output_text_length


def bar(list_value: list[float], array_value: np.ndarray) -> OrderedDict:
    """
    Sample docstring.

    Args:
        list_value (list[float]): _description_
        array_value (np.ndarray): _description_

    Returns:
        OrderedDict: _description_
    """
    assert len(list_value) == array_value.size, "size mismatch"

    output_multiply = np.multiply(np.array(list_value), array_value)

    dict_multiply = OrderedDict()
    for key, value in zip(list_value, output_multiply):
        dict_multiply[key] = value

    return dict_multiply
