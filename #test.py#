from ctypes import *

lib = cdll.LoadLibrary("/home/daniel/Code/D_projects/dag_algorithms/main.so")

func = lib.your_name

func.argtypes = [c_char_p]
func.restype = c_char_p

r = func(b"daniel")
print(r)
