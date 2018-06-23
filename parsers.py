import collections
import typing


class NamedFieldsParser(object):
    """A simple parser that makes all fields accessible by column"""
    def __init__(self, fields: typing.Collection, delimiter='\t'):
        self.fields = fields  # All parsers need this property

        # Internal use
        self._delimiter = delimiter
        self._field_namer = collections.namedtuple('Row', fields)

    def __call__(self, row: str) -> collections.namedtuple:
        row = row.strip().split(self._delimiter)
        return self._field_namer(*row)
