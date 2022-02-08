import sys
import numpy as np
import matplotlib.pyplot as plt
from read_wann_calc_txt import ReadBands, ReadOrbWeight, ReadBerry, ReadOAM
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
	norb = 11

	print("---> Plotting band structure")
	epsk = ReadBands(fname)
	# fout = output_file(out_prefix,"bands")
	# Plot_bandstructure(epsk,klabel)

	print("---> Plotting orbital weight: dz2")
	orb_weight = ReadOrbWeight(fname,norb)
	iorb = 0
	fout = output_file(out_prefix,"weight_dz2")
	Plot_bandstructure_orbweight(epsk,orb_weight,klabel,iorb=iorb,fout=fout)

	print("---> Plotting orbital weight: dx2y2")
	iorb = 3
	fout = output_file(out_prefix,"weight_dx2y2")
	Plot_bandstructure_orbweight(epsk,orb_weight,klabel,iorb=iorb,fout=fout)

	print("---> Plotting orbital weight: dxy")
	iorb = 4
	fout = output_file(out_prefix,"weight_xy")
	Plot_bandstructure_orbweight(epsk,orb_weight,klabel,iorb=iorb,fout=fout)

	print("---> Plotting Berry curvature")
	berry = ReadBerry(fname)
	fout = output_file(out_prefix,"berry_z")
	Plot_bandstructure_berry(epsk,berry,klabel,Bmax=40.0,fout=fout)

	print("---> Plotting OAM (Lz)")
	oam = ReadOAM(fname)
	fout = output_file(out_prefix,"oam_z")
	Plot_bandstructure_oam(epsk,oam,klabel,Lmax=4.0,fout=fout)
#----------------------------------------------------------------------
if __name__ == '__main__':
	main(sys.argv[1:])