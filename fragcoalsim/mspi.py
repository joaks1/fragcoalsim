#! /usr/bin/env python

import os
import glob
import ctypes

_mspi_path_list = glob.glob(os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        "libmspi",
        "libmspi*.so"))
if not _mspi_path_list:
    raise Exception("Could not find bundled libmspi shared library")
_mspi_path = _mspi_path_list[0]
_mspi = ctypes.CDLL(_mspi_path)
_ptr_c_char = ctypes.POINTER(ctypes.c_char)
_ptr_c_double = ctypes.POINTER(ctypes.c_double)
_ptr_ptr_c_char = ctypes.POINTER(_ptr_c_char)
_mspi.run_sims.argtypes = (ctypes.c_int, _ptr_ptr_c_char, _ptr_c_double, _ptr_c_double, _ptr_c_double, ctypes.c_bool)

def mspi_run_sims(args):
    global _mspi
    global _ptr_c_char
    argc = len(args)
    argv = (_ptr_c_char * (argc + 1))()
    for i, arg in enumerate(args):
        enc_arg = arg.encode("utf-8")
        argv[i] = ctypes.create_string_buffer(enc_arg)
    nreps = int(args[2])

    pis = (ctypes.c_double * nreps)()
    pis_between = (ctypes.c_double * nreps)()
    pis_within = (ctypes.c_double * nreps)()
    _mspi.run_sims(argc, argv, pis, pis_between, pis_within, False)
    
    py_pis = [float(x) for x in pis]
    py_pis_between = [float(x) for x in pis_between]
    py_pis_within = [float(x) for x in pis_within]
    return py_pis, py_pis_between, py_pis_within
