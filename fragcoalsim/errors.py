#! /usr/bin/env python

class FragcoalsimError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class ArgumentError(FragcoalsimError):
    def __init__(self, *args, **kwargs):
        FragcoalsimError.__init__(self, *args, **kwargs)

class TempFSError(FragcoalsimError):
    def __init__(self, *args, **kwargs):
        FragcoalsimError.__init__(self, *args, **kwargs)
