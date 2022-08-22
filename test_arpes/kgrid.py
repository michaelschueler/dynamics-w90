import numpy as np
#----------------------------------------------------------------------
def recip_lattice(real_lat):
	recip_lat = np.zeros([3,3])
	recip_lat[1-1, 1-1] = real_lat[2-1, 2-1]*real_lat[3-1, 3-1] - real_lat[3-1, 2-1]*real_lat[2-1, 3-1]
	recip_lat[1-1, 2-1] = real_lat[2-1, 3-1]*real_lat[3-1, 1-1] - real_lat[3-1, 3-1]*real_lat[2-1, 1-1]
	recip_lat[1-1, 3-1] = real_lat[2-1, 1-1]*real_lat[3-1, 2-1] - real_lat[3-1, 1-1]*real_lat[2-1, 2-1]
	recip_lat[2-1, 1-1] = real_lat[3-1, 2-1]*real_lat[1-1, 3-1] - real_lat[1-1, 2-1]*real_lat[3-1, 3-1]
	recip_lat[2-1, 2-1] = real_lat[3-1, 3-1]*real_lat[1-1, 1-1] - real_lat[1-1, 3-1]*real_lat[3-1, 1-1]
	recip_lat[2-1, 3-1] = real_lat[3-1, 1-1]*real_lat[1-1, 2-1] - real_lat[1-1, 1-1]*real_lat[3-1, 2-1]
	recip_lat[3-1, 1-1] = real_lat[1-1, 2-1]*real_lat[2-1, 3-1] - real_lat[2-1, 2-1]*real_lat[1-1, 3-1]
	recip_lat[3-1, 2-1] = real_lat[1-1, 3-1]*real_lat[2-1, 1-1] - real_lat[2-1, 3-1]*real_lat[1-1, 1-1]
	recip_lat[3-1, 3-1] = real_lat[1-1, 1-1]*real_lat[2-1, 2-1] - real_lat[2-1, 1-1]*real_lat[1-1, 2-1]

	vol = real_lat[0,0]*recip_lat[0,0] 
	vol += real_lat[0,1]*recip_lat[0,1]
	vol += real_lat[0,2]*recip_lat[0,2]

	recip_lat = 2*np.pi/vol * recip_lat

	return recip_lat
#----------------------------------------------------------------------
def GenKpath(points,nseg,file_kpts=""):
	kcoords = np.array([points[0]])

	for j in range(1,nseg):
		x = j/nseg
		kp = (1.0-x)*np.array(points[0]) + x*np.array(points[1])
		kcoords = np.append(kcoords,[kp],axis=0)

	for i in range(1,len(points)-1):
		for j in range(nseg):
			x = j/nseg
			kp = (1.0-x)*np.array(points[i]) + x*np.array(points[i+1])
			kcoords = np.append(kcoords,[kp],axis=0)

	kcoords = np.append(kcoords,[points[-1]],axis=0)

	if len(file_kpts) > 0:
		np.savetxt(file_kpts,np.c_[kcoords])

	return kcoords
#----------------------------------------------------------------------
def rtpairs(r, n):

    for i in range(len(r)):
       for j in range(n[i]):    
        yield r[i], j*(2 * np.pi / n[i])
#----------------------------------------------------------------------
def GenKdisk(kpt,krad,Nr,fac,file_kpts):
	# rs = krad * np.sqrt(np.random.random(size=Npts))
	# ts = 2.0 * np.pi * np.random.random(size=Npts)

	rs = krad * np.linspace(0.0,1.0,Nr)
	Ts = np.array(np.floor(fac * rs/krad), dtype=np.int64) + 1

	xs = np.array([])
	ys = np.array([])

	for r, t in rtpairs(rs, Ts):
		xs = np.append(xs, r * np.cos(t) + kpt[0])
		ys = np.append(ys, r * np.sin(t) + kpt[1])

	np.savetxt(file_kpts,np.c_[xs,ys])
#----------------------------------------------------------------------
def GenKgrid(Nu,Nv,file_kpts):
	us = np.linspace(-0.5,0.5,Nu)
	vs = np.linspace(-0.5,0.5,Nv)

	f = open(file_kpts,"w")
	for u in us:
		for v in vs:
			f.write(str(u) + ' ' + str(v) + '\n')
	f.close()
#----------------------------------------------------------------------
def GenKgrid_point(Nu,Nv,u0,v0,rad,file_kpts):
	us = np.linspace(-0.5,0.5,Nu)
	vs = np.linspace(-0.5,0.5,Nv)

	f = open(file_kpts,"w")
	for u in us:
		for v in vs:
			if (u-u0)**2 + (v-v0)**2 < rad**2:
				f.write(str(u) + ' ' + str(v) + '\n')
	f.close()
#----------------------------------------------------------------------
def GenKgrid_square(Nu,Nv,umin,umax,vmin,vmax,file_kpts):
	us = np.linspace(umin, umax, Nu)
	vs = np.linspace(vmin, vmax, Nv)

	f = open(file_kpts,"w")
	for u in us:
		for v in vs:
			f.write(str(u) + ' ' + str(v) + '\n')
	f.close()
#----------------------------------------------------------------------
