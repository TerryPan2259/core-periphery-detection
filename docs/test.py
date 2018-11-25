def division(divident, divisor):
    """
    Division function

    This is an example of function documentation.
    It illustrates how to document parameters, return values
    and their types, and also the exception that a function
    or a module may raise under certain conditions.

    :param divident: operation divident
    :type divident: float
    :param divisor: operation divisor
    :type divisor: float
    :return: division result
    :rtype: float
    :raises ZeroDivisionError: when divisor = 0

    .. note:: This function can accept :class:`int` parameters too.

    .. warning:: ``divisor=0`` will cause :exc:`ZeroDivisionError` exception!

    Example::

        result = division(a, b)
    """
    return divident / divisor
