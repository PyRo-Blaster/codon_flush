import os
import sys

print(os.getcwd())
print(sys.path[0])
print(sys.argv[0])
print(os.path.split(os.path.realpath(__file__))[0])

from util import codon_use
