import sys
import numpy as np
import matplotlib.pyplot as plt
from read_wann_calc import ReadBands, ReadOrbWeight, ReadBerry, ReadOAM
from PlotBands import Plot_bandstructure, Plot_bandstructure_orbweight, \
	Plot_bandstructure_berry, Plot_bandstructure_oam
#----------------------------------------------------------------------
def output_file(prefix,tag):
	if len(prefix) > 0:
		return prefix + "_{}.pdf".format(tag)
	else:
		return ""
#----------------------------------------------------------------------
def main(argv):
	if len(argv) > 0:
		fname = argv[0]
	else:
		print("No input. Exiting ...")
		exit()

	out_prefix = ""
	if len(argv) > 1:
		out_prefix = argv[1]

	klabel = ["M'", "K'", r"$\Gamma$", "K", "M"]

	print("---> Plotting band structure")
	epsk = ReadBands(fname)
	fout = output_file(out_prefix,"bands")
	Plot_bandstructure(epsk,klabel)

	print("---> Plotting orbital weight: dz2")
	orb_weight = ReadOrbWeight(fname)
	iorb = 0
	Plot_bandstructure_orbweight(epsk,orb_weight,klabel,iorb=iorb)

	print("---> Plotting orbital weight: dx2y2")
	orb_weight = ReadOrbWeight(fname)
	iorb = 3
	Plot_bandstructure_orbweight(epsk,orb_weight,klabel,iorb=iorb)

	print("---> Plotting orbital weight: dxy")
	orb_weight = ReadOrbWeight(fname)
	iorb = 4
	Plot_bandstructure_orbweight(epsk,orb_weight,klabel,iorb=iorb)

	print("---> Plotting Berry curvature")
	berry = ReadBerry(fname)
	Plot_bandstructure_berry(epsk,berry,klabel,Bmax=40.0)

	print("---> Plotting OAM (Lz)")
	oam = ReadOAM(fname)
	Plot_bandstructure_oam(epsk,oam,klabel,Lmax=4.0)
#----------------------------------------------------------------------
if __name__ == '__main__':
	main(sys.argv[1:])