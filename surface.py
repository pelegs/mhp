from pylab import *
from prody import *
import argparse
import mhp

parser = argparse.ArgumentParser (description="Calculates MHP for protein surface")
parser.add_argument ('-i','--input', help='Input protein name (without extension)', required=True)
parser.add_argument ('-s','--subfolder', help='Input subfolder name', required=False, default='')
parser.parse_args()
args = vars (parser.parse_args())

pdb = parsePDB(args['subfolder'] + '/' + args['input'] + '.pdb')
psf = parsePSF(args['subfolder'] + '/' + args['input'] + '.psf')
for i in range(25):
    print pdb[i].getCoords(), pdb[i].getElement(), mhp.vdw_radii[pdb[i].getElement()], psf[i].getType(), mhp.F_val[psf[i].getType()]
