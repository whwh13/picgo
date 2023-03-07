from numpy import *
from mayavi import mlab

#SPECIFY THE TRAJ FILE AND THE COLOR FILE
#If no color file is specified, the default coloring will be used	
traj_filename = ''
color_filename = ''

#Open the traj filex
traj_file=open(traj_filename)

#Do we have color?
color=False
if(color_filename != ''):
	color_file = open(color_filename)
	color = True

class FilamentSnapshot:
	def __init__(self):
		self.id=None
		self.type=None
		self.length=None
		self.delta_left=None	
		self.delta_right=None	
		self.coords={}
		self.colors={}
		self.connections=None

class LinkerSnapshot:
	def __init__(self):
		self.id=None
		self.type=None
		self.coords={}
		self.colors={}
		self.connections=None

class MotorSnapshot:
	def __init__(self):
		self.id=None
		self.type=None
		self.coords={}
		self.colors={}
		self.connections=None

class BrancherSnapshot:
	def __init__(self):
		self.id=None
		self.type=None;
		self.coords={}

class BubbleSnapshot:
	def __init__(self):
		self.id=None
		self.type=None;
		self.coords={}

class TrajSnapshot:
	def __init__(self):
		self.step=None
		self.time=None
		self.n_filaments=None
		self.n_linkers=None
		self.n_motors=None
		self.n_branchers=None
		self.n_bubbles=0
		self.filaments={}
		self.linkers={}
		self.motors={}
		self.branchers={}
		self.bubbles={}

TrajSnapshotList=[]
first_snapshot_line=True
first_line=True
reading_filament=False
reading_linker=False
reading_motor=False
reading_brancher=False
reading_bubble=False
n_beads_filament=0
n_beads_linker=0
n_beads_motor=0
n_beads_brancher=0;
n_beads_bubble=0;

line_number = 0;

#read color file
if(color):
	color_lines = [line.strip() for line in color_file]

for line in traj_file:
	line = line.strip()

	if(color):
		color_line = color_lines[line_number]

	if(first_snapshot_line):
		F=TrajSnapshot()

		if len(line.split()) == 7:
			F.step, F.time, F.n_filaments, F.n_linkers, \
			F.n_motors, F.n_branchers, F.n_bubbles = map(double,line.split())
		else:
			F.step, F.time, F.n_filaments, F.n_linkers, \
			F.n_motors, F.n_branchers = map(double,line.split())

		F.n_filaments = int(F.n_filaments)
		F.n_linkers   = int(F.n_linkers)
		F.n_motors    = int(F.n_motors)
		F.n_branchers = int(F.n_branchers)
		F.n_bubbles   = int(F.n_bubbles)
		first_snapshot_line=False
		line_number+=1
		continue

	if(len(line)==0):
		first_snapshot_line=True
		first_filament_line=True
		TrajSnapshotList.append(F)
		assert F.n_filaments == len(F.filaments)
		assert F.n_linkers   == len(F.linkers)
		assert F.n_motors    == len(F.motors)
		assert F.n_branchers == len(F.branchers)
		assert F.n_bubbles   == len(F.bubbles)
		n_beads_filament=0
		n_beads_linker=0
		n_beads_motor=0
		n_beads_brancher=0;
		n_beads_bubble=0;
		line_number+=1
		continue
			
	if(first_line):
		line_split = line.split()

		if(line_split[0] == "FILAMENT") or (line_split[0] == "F"):
			FS=FilamentSnapshot()
			
			if(len(line_split[1:]) == 5):
				FS.id, FS.type, FS.length, FS.delta_left, FS.delta_right = map(int,line_split[1:])
			else:
				FS.id, FS.length, FS.delta_left, FS.delta_right = map(int,line_split[1:])

			first_line=False
			F.filaments[FS.id]=FS
			FS.connections=vstack(
				[arange(n_beads_filament, n_beads_filament + FS.length - 1.5), 
				arange(n_beads_filament + 1, n_beads_filament + FS.length - .5)]).T
			n_beads_filament+=FS.length
			reading_filament=True
			line_number+=1
			continue
		if(line_split[0] == "LINKER") or (line_split[0] == "L"):
			LS = LinkerSnapshot()
			LS.id, LS.type = map(int, line_split[1:])
			first_line = False
			F.linkers[LS.id] = LS
			LS.connections=vstack([arange(n_beads_linker, n_beads_linker + 0.5), 
				arange(n_beads_linker + 1, n_beads_linker + 1.5)]).T
			n_beads_linker+=2
			reading_linker=True
			line_number+=1
			continue
		if(line_split[0] == "MOTOR") or (line_split[0] == "M"):
			MS = MotorSnapshot()
			MS.id, MS.type = map(int, line_split[1:])
			first_line = False
			F.motors[MS.id] = MS
			MS.connections=vstack([arange(n_beads_motor, n_beads_motor + 0.5), 
				arange(n_beads_motor + 1, n_beads_motor + 1.5)]).T
			n_beads_motor+=2
			reading_motor=True
			line_number+=1
			continue
		if(line_split[0] == "BRANCHER") or (line_split[0] == "B"):
			BS = BrancherSnapshot()
			BS.id, BS.type = map(int, line_split[1:])
			first_line = False
			F.branchers[BS.id] = BS
			n_beads_brancher+=1
			reading_brancher=True
			line_number+=1
			continue
		if(line_split[0] == "BUBBLE"):
			US = BubbleSnapshot()
			US.id, US.type = map(int, line_split[1:])
			first_line = False
			F.bubbles[US.id] = US
			n_beads_bubble+=1
			reading_bubble=True
			line_number+=1
			continue

	if(reading_filament):
		FS.coords=array(line.split(),'f')
		if(color):
			FS.colors=array(color_line.split(), 'f')
		N=len(FS.coords)
		FS.coords=FS.coords.reshape(N/3,3).T
		first_line=True
		reading_filament=False

	if(reading_linker):
		LS.coords=array(line.split(),'f')
		if(color):
			LS.colors=array(color_line.split(), 'f')
		N=len(LS.coords)
		LS.coords=LS.coords.reshape(N/3,3).T
		first_line=True
		reading_linker=False

	if(reading_motor):
		MS.coords=array(line.split(),'f')
		if(color):
			MS.colors=array(color_line.split(), 'f')
		N=len(MS.coords)
		MS.coords=MS.coords.reshape(N/3,3).T
		first_line=True
		reading_motor=False

	if(reading_brancher):
		BS.coords=array(line.split(),'f')
		N=len(BS.coords)
		BS.coords=BS.coords.reshape(N/3,3).T
		first_line=True
		reading_brancher=False

	if(reading_bubble):
		US.coords=array(line.split(),'f')
		N=len(US.coords)
		US.coords=US.coords.reshape(N/3,3).T
		first_line=True
		reading_bubble=False

	line_number+=1


#Code to implicitly plot an algebraic surface, by:
#http://indranilsinharoy.com/2014/03/02/plotting-algebraic-surfaces-using-mayavi/
def implicit_plot(expr, ext_grid, fig_handle=None, Nx=101, Ny=101, Nz=101,
                 col_isurf=(50/255, 199/255, 152/255), col_osurf=(240/255,36/255,87/255),
                 opa_val=0.8, opaque=True, ori_axis=True, **kwargs):

    if fig_handle==None:  # create a new figure
        fig = mlab.figure(1,bgcolor=(0.97, 0.97, 0.97), fgcolor=(0, 0, 0), size=(800, 800))
    else:
        fig = fig_handle

    xl, xr, yl, yr, zl, zr = ext_grid
    x, y, z = np.mgrid[xl:xr:eval('{}j'.format(Nx)),
                       yl:yr:eval('{}j'.format(Ny)),
                       zl:zr:eval('{}j'.format(Nz))]
    scalars = eval(expr)
    src = mlab.pipeline.scalar_field(x, y, z, scalars)
    if opaque:
        delta = 1.e-5
        opa_val=1.0
    else:
        delta = 0.0

    cont1 = mlab.pipeline.iso_surface(src, color=col_isurf, contours=[0-delta],
                                      transparent=False, opacity=opa_val)
    cont1.compute_normals = False # for some reasons, setting this to true actually cause
                                  # more unevenness on the surface, instead of more smooth

    if opaque: # the outer surface is specular, the inner surface is not
        cont2 = mlab.pipeline.iso_surface(src, color=col_osurf, contours=[0+delta],
                                          transparent=False, opacity=opa_val)
        cont2.compute_normals = False
        cont1.actor.property.backface_culling = True
        cont2.actor.property.frontface_culling = True
        cont2.actor.property.specular = 0.2 #0.4 #0.8
        cont2.actor.property.specular_power = 15.0

    else:  # make the surface (the only surface) specular
        cont1.actor.property.specular = 0.2 #0.4 #0.8
        cont1.actor.property.specular_power = 15.0
 
    # Scene lights (4 lights are used)
    engine = mlab.get_engine()
    scene = engine.current_scene
    cam_light_azimuth = [78, -57, 0, 0]
    cam_light_elevation = [8, 8, 40, -60]
    cam_light_intensity = [0.72, 0.48, 0.60, 0.20]

    for i in range(4):
        camlight = scene.scene.light_manager.lights[i]
        camlight.activate = True
        camlight.azimuth = cam_light_azimuth[i]
        camlight.elevation = cam_light_elevation[i]
        camlight.intensity = cam_light_intensity[i]

    # axis through the origin
    if ori_axis:
        len_caxis = int(1.05*np.max(np.abs(np.array(ext_grid))))
        caxis = mlab.points3d(0.0, 0.0, 0.0, len_caxis, mode='axes',color=(0.15,0.15,0.15),
                              line_width=1.0, scale_factor=1.,opacity=1.0)
        caxis.actor.property.lighting = False

    # if no figure is passed, the function will create a figure.
    if fig_handle==None:
        # Setting camera
        cam = fig.scene.camera
        cam.elevation(-20)
        cam.zoom(1.0) # zoom should always be in the end.
        mlab.show()


@mlab.show
def show_snapshot(snapshot_number=-1):

	#if were saving the Snapshots
	saving = False
	saveFile = ''

	#PARAMETERS TO SET FOR VISUAL
	#for color scaling
	MAXVAL = 0
	MINVAL = 0
	SCALETITLE = ''
	COLORMAP = ''

	#tracking plus and minus ends
	TRACKENDS = False

	#color by radial angle
	COLORBYANGLE = False

	#grid size
	GRIDSIZEMAXX = 0.0
	GRIDSIZEMINX = 0.0

	GRIDSIZEMAXY = 0.0
	GRIDSIZEMINY = 0.0

	GRIDSIZEMAXZ = 0.0
	GRIDSIZEMINZ = 0.0

	#boundary type, CUBIC or SPHERICAL
	BOUNDARYTYPE = "CUBIC"

	#diameter if spherically bound
	DIAMETER = 0.0

	#default color, in RGB
	DBEADCOLOR    = (1.0,1.0,1.0) 
	DFILCOLOR     = (1.0,0.0,0.0)
	DLINKERCOLOR  = (0.0,1.0,0.1)
	DMOTORCOLOR   = (0.0,0.2,1.0)
	DBUBBLECOLOR  = (0.2,0.7,0.5)
	PLUSENDCOLOR  = (0.0,0.0,0.0)
	MINUSENDCOLOR = (1.0,1.0,1.0)

	local_snapshot=TrajSnapshotList[snapshot_number]
	figw = mlab.figure(1, size=(1000, 1000), bgcolor=(1.0,1.0,1.0))
	mlab.clf()

	#create grid
	if BOUNDARYTYPE == "CUBIC" :

		x = [GRIDSIZEMINX,GRIDSIZEMINX,GRIDSIZEMAXX]
		y = [GRIDSIZEMAXY,GRIDSIZEMINY,GRIDSIZEMINY]
		z = [GRIDSIZEMINZ,GRIDSIZEMAXZ,GRIDSIZEMINZ]
		pts = mlab.points3d(x,y,z, scale_mode='none', scale_factor=0.2)
		
		outline=mlab.pipeline.outline(pts, line_width=0.6, color=(0.0,0.0,0.0))

	elif BOUNDARYTYPE == "SPHERICAL" :

		implicit_plot('(x - {a})**2 + (y - {b})**2 + (z - {c})**2 - {r}**2'.format(a=GRIDSIZEMAXX/2, b=GRIDSIZEMAXY/2, c=GRIDSIZEMAXZ/2, r=DIAMETER/2), 
					  (0, GRIDSIZEMAXX, 0, GRIDSIZEMAXY, 0, GRIDSIZEMAXZ),
              		  fig_handle=figw, Nx=64, Ny=64, Nz=64, col_isurf=(0.67, 0.77, 0.93),
             		  opaque=False, opa_val=0.3, ori_axis=False)

	elif BOUNDARYTYPE == "ELLIPSOID" :

		implicit_plot('((x - {i})/ {a})**2 + ((y - {j}) / {b})**2 + ((z - {k})/ {c})**2 - 1'.format(a=GRIDSIZEMAXX/2, b=GRIDSIZEMAXY/2, c=GRIDSIZEMAXZ/2, 
																									i=GRIDSIZEMAXX/2, j=GRIDSIZEMAXY/2, k=GRIDSIZEMAXZ/2), 
					  (0 - 200, GRIDSIZEMAXX + 200, 0 - 200, GRIDSIZEMAXY + 200, 0 - 200, GRIDSIZEMAXZ + 200),
              		  fig_handle=figw, Nx=64, Ny=64, Nz=64, col_isurf=(0.67, 0.77, 0.93),
             		  opaque=False, opa_val=0.3, ori_axis=False)
	#display time
	time = 'Time = ' + str(local_snapshot.time) + "s"
	mlab.text(0.6, 0.9, time, color=(0.0,0.0,0.0))

	#DISPLAYING RANDOM POINTS FOR MONOMERS
	#can add as many types of monomers as needed
	n1_monomers = 1000

	n1_x = []
	n1_y = []
	n1_z = []

	if(BOUNDARYTYPE == "CUBIC"):

		for i in xrange(0, n1_monomers):

			n1_x.append(random.uniform(0,GRIDSIZEMAXX))
			n1_y.append(random.uniform(0,GRIDSIZEMAXY))
			n1_z.append(random.uniform(0,GRIDSIZEMAXZ))

		mlab.points3d(n1_x,n1_y,n1_z, scale_factor=12.0, color=DFILCOLOR)

	elif BOUNDARYTYPE == "SPHERICAL" :

		for i in xrange(0, n1_monomers):

			r = random.uniform(0,DIAMETER/2)
			l = random.uniform(0,2 * math.pi)
			h = random.uniform(-math.pi/2, math.pi/2)

			n1_x.append(r * math.cos(l) * math.cos(h) + GRIDSIZEMAXX/2)
			n1_y.append(r * math.sin(h) + GRIDSIZEMAXY/2)
			n1_z.append(r * math.sin(l) * math.cos(h) + GRIDSIZEMAXZ/2)

		mlab.points3d(n1_x,n1_y,n1_z, scale_factor=12.0, color=DFILCOLOR)

	#DISPLAYING FILAMENTS
	if(len(local_snapshot.filaments) != 0):

		plusends =[[], [], []]
		minusends=[[], [], []]

		#can add as many filaments as desired
		x=[]
		c=[]
		connections=[]

		for i in 0, 1, 2:
			q=[]
			for fid in sorted(local_snapshot.filaments.keys()):
				q.append(local_snapshot.filaments[fid].coords[i])

			x.append(hstack(q))

		#plus ends and minus ends if tracked
		if(TRACKENDS):
			for i in 0, 1, 2:
				for fid in sorted(local_snapshot.filaments.keys()):
					coords = local_snapshot.filaments[fid].coords

					plusends[i].append(coords[i][len(coords[i]) - 1])
					minusends[i].append(coords[i][0])

		for fid in sorted(local_snapshot.filaments.keys()):
			for color in local_snapshot.filaments[fid].colors:
				c.append(color)

		if (COLORBYANGLE): 
			for filamentID, FS in local_snapshot.filaments.iteritems():

				#loop through cylinders
				for i in xrange(1, len(FS.coords[0]), 1):

					b1x = FS.coords[0][i-1]
					b1y = FS.coords[1][i-1]
					b1z = FS.coords[2][i-1]

					b2x = FS.coords[0][i]
					b2y = FS.coords[1][i]
					b2z = FS.coords[2][i]

					b1 = np.array([b1x, b1y, b1z])
					b2 = np.array([b2x, b2y, b2z])

					#calculate angle from z plane
					dist = sqrt((b2 - b1).dot(b2 - b1))
					direction = (b2 - b1) / dist

					angle = math.acos(direction.dot(np.array([0,0,1]))) * 180 / math.pi - 90

					c.append(-angle)

					if (i == len(FS.coords[0]) - 1): 
						c.append(angle)

		for fid in sorted(local_snapshot.filaments.keys()):
			connections.append(local_snapshot.filaments[fid].connections)

		connections = vstack(connections)

		# Create the points
		if(len(c) != 0):
			src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2], c)
		else:
			src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2])

		if(TRACKENDS):
			plusendsrc = mlab.pipeline.scalar_scatter(plusends[0], plusends[1], plusends[2])
			gsphere=mlab.pipeline.glyph(plusendsrc, mode="sphere", resolution=24, 
								        scale_mode='scalar', scale_factor=20.0, color=(PLUSENDCOLOR))

			minusendsrc = mlab.pipeline.scalar_scatter(minusends[0], minusends[1], minusends[2])
			gsphere=mlab.pipeline.glyph(minusendsrc, mode="sphere", resolution=24, 
								        scale_mode='scalar', scale_factor=20.0, color=(MINUSENDCOLOR))

		gsphere=mlab.pipeline.glyph(src, mode="sphere", resolution=24, 
								     scale_mode='scalar', scale_factor=3.0, color=(DBEADCOLOR))

		# Connect them
		src.mlab_source.dataset.lines = connections

		# Finally, display the set of lines
		tube=mlab.pipeline.tube(src, tube_radius=5)
		tube.filter.number_of_sides=12

		if(len(c) != 0):
			surface = mlab.pipeline.surface(tube, colormap=COLORMAP,
										    vmax = MAXVAL, vmin = MINVAL)
			cb = mlab.colorbar(object=surface, orientation='vertical', 
						  	   title=SCALETITLE, label_fmt='%.0f', nb_labels=5)

			cb.title_text_property.color = (0.0,0.0,0.0)
			cb.label_text_property.color = (1.0,1.0,1.0)

		else:
			surface = mlab.pipeline.surface(tube, color=DFILCOLOR)

	#DISPLAYING LINKERS
	if(len(local_snapshot.linkers) != 0):
		x=[]
		c=[]
		connections=[]

		#Add random diffusing linkers
		n_diffusing_linkers = 0
		len_linker = 35.0
		linker_type = 0
		n_beads_linker = 2 * len(local_snapshot.linkers) 

		for i in xrange(n_beads_linker, n_beads_linker + 2 * n_diffusing_linkers, 2) :

			LS = LinkerSnapshot()
			LS.id = random.randint(10E6,10E8)
			LS.type = linker_type
			local_snapshot.linkers[LS.id] = LS
			LS.connections=vstack([arange(i, i + 0.5), 
								   arange(i + 1, i + 1.5)]).T

			#first random coord
			if BOUNDARYTYPE == "CUBIC":

				coord1x = random.uniform(0,GRIDSIZEMAXX - len_linker)
				coord1y = random.uniform(0,GRIDSIZEMAXY - len_linker)
				coord1z = random.uniform(0,GRIDSIZEMAXZ - len_linker)

			elif BOUNDARYTYPE == "SPHERICAL":

				r = random.uniform(0,DIAMETER/2)
				l = random.uniform(0,2 * math.pi)
				h = random.uniform(-math.pi/2, math.pi/2)

				coord1x = r * math.cos(l) * math.cos(h) + GRIDSIZEMAXX/2
				coord1y = r * math.sin(h) + GRIDSIZEMAXY/2
				coord1z = r * math.sin(l) * math.cos(h) + GRIDSIZEMAXZ/2

			#random projection
			direction = [random.uniform(-1,1),
						 random.uniform(-1,1),
						 random.uniform(-1,1)]

			dist = sqrt(pow(direction[0],2) + 
						pow(direction[1],2) + 
						pow(direction[2],2))

			normDir = [float(i)/dist for i in direction]

			coord2x = coord1x + normDir[0] * len_linker
			coord2y = coord1y + normDir[1] * len_linker
			coord2z = coord1z + normDir[2] * len_linker

			LS.coords=array([coord1x, coord1y, coord1z,
					   		 coord2x, coord2y, coord2z])

			N=len(LS.coords)
			LS.coords=LS.coords.reshape(N/3,3).T

		for i in 0, 1, 2:
			q=[]
			for lid in sorted(local_snapshot.linkers.keys()):
				q.append(local_snapshot.linkers[lid].coords[i])

			x.append(hstack(q))

		for lid in sorted(local_snapshot.linkers.keys()):
			for color in local_snapshot.linkers[lid].colors:
				c.append(color)

		for lid in sorted(local_snapshot.linkers.keys()):
			connections.append(local_snapshot.linkers[lid].connections)

		connections = vstack(connections)

		# Create the points
		if(len(c) != 0):
			src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2], c)
		else:
			src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2])

		gsphere=mlab.pipeline.glyph(src, mode="sphere", resolution=24, 
								    scale_mode='none', color=(DBEADCOLOR))

		# Connect them
		src.mlab_source.dataset.lines = connections

		# Finally, display the set of lines
		tube=mlab.pipeline.tube(src, tube_radius=3)
		tube.filter.number_of_sides=12

		if(len(c) != 0):
			surface = mlab.pipeline.surface(tube, colormap=COLORMAP, 
											vmax = MAXVAL, vmin = MINVAL)
		else:
			surface = mlab.pipeline.surface(tube, color=DLINKERCOLOR)

	#DISPLAYING MOTORS
	if(len(local_snapshot.motors) != 0):
		x=[]
		c=[]
		connections=[]

		#Add random diffusing motors
		n_diffusing_motors = 0
		len_motor = 200.0
		motor_type = 0
		n_beads_motor = 2 * len(local_snapshot.motors) 

		for i in xrange(n_beads_motor, n_beads_motor + 2 * n_diffusing_motors, 2) :

			MS = MotorSnapshot()
			MS.id = random.randint(10E6,10E8)
			MS.type = motor_type
			local_snapshot.motors[MS.id] = MS
			MS.connections=vstack([arange(i, i + 0.5), 
								   arange(i + 1, i + 1.5)]).T

			#first random coord
			if BOUNDARYTYPE == "CUBIC":

				coord1x = random.uniform(0,GRIDSIZEMAXX - len_motor)
				coord1y = random.uniform(0,GRIDSIZEMAXY - len_motor)
				coord1z = random.uniform(0,GRIDSIZEMAXZ - len_motor)

			elif BOUNDARYTYPE == "SPHERICAL":

				r = random.uniform(0,DIAMETER/2)
				l = random.uniform(0,2 * math.pi)
				h = random.uniform(-math.pi/2, math.pi/2)

				coord1x = r * math.cos(l) * math.cos(h) + GRIDSIZEMAXX/2
				coord1y = r * math.sin(h) + GRIDSIZEMAXY/2
				coord1z = r * math.sin(l) * math.cos(h) + GRIDSIZEMAXZ/2

			#random projection
			direction = [random.uniform(-1,1),
						 random.uniform(-1,1),
						 random.uniform(-1,1)]

			dist = sqrt(pow(direction[0],2) + 
						pow(direction[1],2) + 
						pow(direction[2],2))

			normDir = [float(i)/dist for i in direction]

			coord2x = coord1x + normDir[0] * len_motor
			coord2y = coord1y + normDir[1] * len_motor
			coord2z = coord1z + normDir[2] * len_motor

			MS.coords=array([coord1x, coord1y, coord1z,
					   		 coord2x, coord2y, coord2z])

			N=len(MS.coords)
			MS.coords=MS.coords.reshape(N/3,3).T

		for i in 0, 1, 2:
			q=[]
			for mid in sorted(local_snapshot.motors.keys()):
				q.append(local_snapshot.motors[mid].coords[i])
			x.append(hstack(q))

		for mid in sorted(local_snapshot.motors.keys()):
			for color in local_snapshot.motors[mid].colors:
				c.append(color)

		for mid in sorted(local_snapshot.motors.keys()):
				connections.append(local_snapshot.motors[mid].connections)

		connections = vstack(connections)

		# Create the points
		if(len(c) != 0):
			src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2], c)
		else:
			src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2])

		gsphere=mlab.pipeline.glyph(src, mode="sphere", resolution=24, 
									scale_mode='none', color=DBEADCOLOR)

		# Connect them
		src.mlab_source.dataset.lines = connections

		# Finally, display the set of lines
		tube=mlab.pipeline.tube(src, tube_radius=4)
		tube.filter.number_of_sides=12

		if(len(c) != 0):
			surface = mlab.pipeline.surface(tube, colormap=COLORMAP, 
											vmax = MAXVAL, vmin = MINVAL)
		else:
			surface = mlab.pipeline.surface(tube, color=DMOTORCOLOR)

	#DISPLAYING BRANCHERS (NO COLOR)
	if(len(local_snapshot.branchers) != 0):
		x=[]

		for i in 0, 1, 2:
			q=[]
			for bid in sorted(local_snapshot.branchers.keys()):
				q.append(local_snapshot.branchers[bid].coords[i])
			x.append(hstack(q))

		# Create the points
		src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2])
		gsphere=mlab.pipeline.glyph(src, mode="sphere", resolution=24, 
									scale_factor=2.0, color=DBEADCOLOR)

	#DISPLAYING BUBBLES
	if(len(local_snapshot.bubbles) != 0):
		x=[]

		for i in 0, 1, 2:
			q=[]
			for uid in sorted(local_snapshot.bubbles.keys()):
				q.append(local_snapshot.bubbles[uid].coords[i])
			x.append(hstack(q))

		# Create the points
		src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2])
		gsphere=mlab.pipeline.glyph(src, mode="sphere", resolution=24, 
									scale_factor=200.0, color=DBUBBLECOLOR)

	if(saving):
		mlab.savefig(filename=saveFile + "Snapshot" + str(snapshot_number) + ".png")

@mlab.animate(delay=10, ui=True)
def anim():
	for i in range(0,len(TrajSnapshotList), 1):
		show_snapshot(i)
		yield
		

