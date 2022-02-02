import sys
import numpy as np
#--------------------------------------------------------------------------------------
def Gen_Sin2_pulse(F0,w0,nc,t0,nt=400):
	Tp = 2*np.pi/w0 * nc
	ts = np.linspace(t0, t0 + Tp,nt)
	f = F0 * np.sin(np.pi*(ts-t0)/Tp)**2 * np.sin(w0*(ts-t0))
	return ts,f
#--------------------------------------------------------------------------------------
def main(argv):
	Ry = 27.211386 # Rydberg constant

	w0 = 1.5 # eV
	F0 = 1.0e-3
	nc = 6
	t0 = 0.0

	pol = np.array([1.0, 0.0, 0.0])

	file_field = "inp/sin2_pulse_F{}_w{}_nc{}_t{}.dat".format(F0,w0,nc,t0)
	ts, pulse = Gen_Sin2_pulse(F0,w0/Ry,nc,t0,nt=400)
	np.savetxt(file_field, np.c_[ts, pol[0]* pulse, pol[1] * pulse, pol[2] * pulse])
#--------------------------------------------------------------------------------------
if __name__ == '__main__':
	main(sys.argv[1:])