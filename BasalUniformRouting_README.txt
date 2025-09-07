
Data for the Greenland-wide subglacial water system with uniform basal meltrate
load with numpy.load('filename.npy')

BasalUniformRouting.npy:
	numpy arra of shape 6 x 67562 x 3
	first axis is the floatation fraction, six rows for six FFs [0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
	second axis is the mesh, 67562 is the number of horizontal mesh vertices, mesh coordinates in file coordinates2d_xy.txt
	third axis contains the three variable fields bua, bvol, and bflux

BasalUniformRoutingRvVw.npy:
	numpy array of shape 6 x 67562 x 2
	first and second axis same as for BasalUniformRouting.npy
	third axis has two slices of dtype=object
	first slice contains list of receiving vertices for each floatation, second slice contains respective vertex weights
	each entry being an np.array of variable length that lists the vertex numbers/weights (hence dtype=object)
