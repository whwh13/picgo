from numpy import *
import pylab as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import scipy.linalg as la
from scipy.stats import linregress
import math

#Various snapshots for objects
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

#A cylinder snapshot
class CylinderSnapshot:
	def __init__(self):
		self.id = None
		self.coords={}

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

#A full trajectory snapshot
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
		self.cylinders={}
		self.linkers={}
		self.motors={}
		self.branchers={}
		self.bubbles={}


#A chemical snapshot
class ChemSnapshot:
	def __init__(self):
		self.step=None
		self.time=None
		self.diffusingSpecies={}
		self.bulkSpecies={}
		self.filamentSpecies={}
		self.plusEndSpecies={}
		self.minusEndSpecies={}
		self.linkerSpecies={}
		self.motorSpecies={}
		self.brancherSpecies={}


#A histogram of values from a snapshot
class HistSnapshot:
	def __init__(self):
		self.step=None
		self.time=None
		self.hist={}

########################################################
#
#			         PARSING FUNCTIONS
#
########################################################


#Read all data from a trajectory file
#returns a list of Snapshots with all data
def readTrajectory(filename):

	print "Reading " + filename + "..."

	#Open the traj file
	traj_file=open(filename)

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

	for line in traj_file:
		line = line.strip()

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
			N=len(FS.coords)
			FS.coords=FS.coords.reshape(N/3,3)

			#create cylinder snapshots
			for i in xrange(0, len(FS.coords) - 1):
				CS = CylinderSnapshot()
				#create unique id using cylinder position and fil id
				CS.id = (FS.id + i) * (FS.id + i + 1) / 2 + FS.id
				CS.coords = [FS.coords[i], FS.coords[i+1]]
				F.cylinders[CS.id] = CS

			first_line=True
			reading_filament=False

		if(reading_linker):
			LS.coords=array(line.split(),'f')
			N=len(LS.coords)
			LS.coords=LS.coords.reshape(N/3,3)
			first_line=True
			reading_linker=False

		if(reading_motor):
			MS.coords=array(line.split(),'f')
			N=len(MS.coords)
			MS.coords=MS.coords.reshape(N/3,3)
			first_line=True
			reading_motor=False

		if(reading_brancher):
			BS.coords=array(line.split(),'f')
			N=len(BS.coords)
			BS.coords=BS.coords.reshape(N/3,3)
			first_line=True
			reading_brancher=False

		if(reading_bubble):
			US.coords=array(line.split(),'f')
			N=len(US.coords)
			US.coords=US.coords.reshape(N/3,3)
			first_line=True
			reading_bubble=False

		line_number+=1

	#return the final Snapshot list
	print str(len(TrajSnapshotList)) + " snapshots read."

	return TrajSnapshotList	

#Read a number of trajectories, return array of Snapshot lists
def readTrajectories(fileNames, fileDirectory='', runMin=0, runMax=10):

	#if file directory specified, auto-generate file names:
	if (fileDirectory != ''):

		fileNames = []

		for trajNum in xrange(runMin, runMax):
			fileNames.append(fileDirectory + '/Run' + str(trajNum) + '/snapshot.traj')

	#now read data
	TrajSnapshotLists = []

	for trajFile in fileNames:
		print "Reading trajectory..."
		TrajSnapshotLists.append(readTrajectory(trajFile))

	return TrajSnapshotLists

#Read all data from a chemical output file
#return a list of chemical snapshots
def readChemistry(fileName):

	print "Reading " + fileName + "..."

	ChemSnapshotList = []

	#Open the traj file
	chem_file=open(fileName)

	first_Snapshot_line=True
	first_line=True
	line_number=0

	for line in chem_file:
		line = line.strip()

		if(first_Snapshot_line):
			C=ChemSnapshot()
			C.step, C.time = map(double,line.split())
			first_Snapshot_line=False
			line_number+=1
			continue

		elif(len(line)==0):
			first_Snapshot_line=True
			ChemSnapshotList.append(C)
			line_number+=1
			continue

		else:
			line_split = line.split()
			
			species_split = line_split[0].split(":")

			species_name = species_split[0]
			species_type = species_split[1]

			if(species_type == "DIFFUSING"):
				C.diffusingSpecies[species_name] = line_split[1]

			if(species_type == "BULK"):
				C.bulkSpecies[species_name] = line_split[1]

			if(species_type == "FILAMENT"):
				C.filamentSpecies[species_name] = line_split[1]

			if(species_type == "PLUSEND"):
				C.plusEndSpecies[species_name] = line_split[1]

			if(species_type == "MINUSEND"):
				C.minusEndSpecies[species_name] = line_split[1]

			if(species_type == "LINKER"):
				C.linkerSpecies[species_name] = line_split[1]

			if(species_type == "MOTOR"):
				C.motorSpecies[species_name] = line_split[1]

			if(species_type == "BRANCHER"):
				C.brancherSpecies[species_name] = line_split[1]

			line_number+=1

	print str(len(ChemSnapshotList)) + " snapshots read."

	return ChemSnapshotList


#Read a number of chem outputs, return array of chem lists
def readChemistry(fileNames):

	ChemSnapshotLists = []

	for chemFile in fileNames:
		print "Reading chemistry output..."
		ChemSnapshotLists.append(readChemistryOutput(chemFile))

	return ChemSnapshotLists

#Read a histogram from a file, return histogram snapshot list
def readHistogram(fileName):

	print "Reading " + fileName + "..."

	HistSnapshotList = []

	#Open the traj file
	hist_file=open(fileName)

	first_Snapshot_line=True
	first_line=True
	line_number=0

	for line in hist_file:
		line = line.strip()

		if(first_Snapshot_line):
			H=HistSnapshot()
			H.step, H.time = map(double,line.split())
			first_Snapshot_line=False
			line_number+=1
			continue

		elif(len(line)==0):
			first_Snapshot_line=True
			HistSnapshotList.append(H)
			line_number+=1
			continue

		else:
			line_split = line.split()
			
			for i in xrange(0, len(line_split), 3):
				binMin, binMax, freq = map(double, line_split[i:i+3])
				H.hist[(binMin, binMax)] = freq

			line_number+=1

	print str(len(HistSnapshotList)) + " snapshots read."

	return HistSnapshotList

#Read a number of histogram outputs, return array of histogram lists
def readHistograms(fileNames):

	HistSnapshotLists = []

	for histFile in fileNames:
		print "Reading histogram output..."
		HistSnapshotLists.append(readHistograms(histFile))

	return HistogramSnapshotLists

#combine two histograms
#this can be called recursively to get a single histogram for set of trajectories
#returns a single histogram snapshot
def combineTwoHistograms(histSnapshot1, histSnapshot2):

	H = HistSnapshot()

	#avg time and step
	H.time = (histSnapshot1.time + histSnapshot2.time) / 2
	H.step = (histSnapshot1.step + histSnapshot2.step) / 2

	#combine hist values
	for bins, freq in histSnapshot1.hist.iteritems():
		H.hist[bins] = freq

	for bins, freq in histSnapshot2.hist.iteritems():
		H.hist[bins] += freq

	return H


#combine an arbitrary length histogram snapshot list
def combineAllHistograms(histSnapshotList):

	#combine first two
	H = combineTwoHistograms(histSnapshotList[0], histSnapshotList[1])

	for i in xrange(2, len(histSnapshotList), 1):
		H = combineTwoHistograms(H, histSnapshotList[i])

	return H

########################################################
#
#
#			     ANALYSIS FUNCTIONS
#
#	contains: 
#		-average filament end to end distance
#		-average filament length
#		-minimum sphere enclosing volume
#		-Nematic order parameter
#	    -contractile velocity
#		-density of entire system
#		
#
########################################################


#Calculate average filament end to end distance at each step
#returns a time and distance pairs
def avgEndToEndDistance(SnapshotList, snapshot=1):

	Snapshot = SnapshotList[snapshot]

	totalDist = 0;

	for filamentID, FS in Snapshot.filaments.iteritems():

		#get first and last bead
		b1 = FS.coords[0]
		b2 = FS.coords[len(FS.coords) - 1]

		#calculate distance
		dist = sqrt((b2 - b1).dot(b2 - b1))
		totalDist += dist

	#average
	totalDist /= len(Snapshot.filaments)
	totalDist *= 0.001 #to um

	#add to list
	return (Snapshot.time, totalDist)


#Calculate average filament end to end distance from a number of trajectories
def avgEndToEndDistances(SnapshotLists, snapshot=1):

	avgDistsList = []

	#get all data
	for SnapshotList in SnapshotLists:

		avgDists = avgEndToEndDistance(SnapshotList, snapshot)
		avgDistsList.append(avgDists)

	#now combine data
	finalTime = 0
	finalAvgDist = 0

	for avgDists in avgDistsList:

		finalTime += avgDists[0] 
		finalAvgDist += avgDists[1]

	#average
	finalTime /= len(SnapshotLists)
	finalAvgDist /= len(SnapshotLists)
	error = numpy.std(avgDistsList)

	return (finalTime, finalAvgDist, error)

#Calculate average filament length at each step
#returns a time and distance pair
def avgLength(SnapshotList, snapshot=1):

	Snapshot = SnapshotList[snapshot]

	totalLength = 0;

	for filamentID, FS in Snapshot.filaments.iteritems():

		filLength = 0;

		for i in xrange(0, len(FS.coords) - 1):

			#get first and last bead
			b1 = FS.coords[i]
			b2 = FS.coords[i+1]

			#calculate distance
			filLength += sqrt((b2 - b1).dot(b2 - b1))
			
		totalLength += filLength

	#average
	totalLength /= len(Snapshot.filaments)

	#add to list
	return (Snapshot.time,totalLength)

#Calculate average lengths from a number of trajectories
def avgLengths(SnapshotLists, snapshot=1):

	avgLengthsList = []

	#get all data
	for SnapshotList in SnapshotLists:

		avgLengths = avgLength(SnapshotList, snapshot)
		avgLengthsList.append(avgLengths)

	#now combine data
	finalTime = 0
	finalAvgLength = 0

	for avgLengths in avgLengthsList:

		finalTime += avgLengths[0]
		finalAvgLength += avgLengths[1]

	#average
	finalTime /= len(SnapshotLists)
	finalAvgLength /= len(SnapshotLists)

	error = numpy.std(avgLengthsList)

	return (finalTime, finalAvgLength, error)

#CREDIT TO: BlenderArtists.org

#minimum covering sphere candidate (mcsc)
#given candidate_points find minimum covering sphere
#find coordinates of point(candidate center) that is equidistant from candidate points
#if a coordinate is negative, point(candidate center) is outside of simplex,
#and the candidate point with the smallest value is removed

def mcsc(points, candidate_indices):

    candidate_count = len(candidate_indices)
    candidate_center = np.array([0.0,0.0,0.0])

    if candidate_count == 1:

        candidate_center = points[candidate_indices[0]]

    if candidate_count == 2:

        q0 = points[candidate_indices[0]]
        q1 = points[candidate_indices[1]]
        #l0 = l1 = 0.5
        #candidate_center then is midpoint formula
        candidate_center = (q0 + q1)/2.0

    if candidate_count == 3:

        q0 = points[candidate_indices[0]]
        q1 = points[candidate_indices[1]]
        q2 = points[candidate_indices[2]]

        a00 = (q0 - q2).dot(q0 - q2)
        a01 = (q0 - q2).dot(q1 - q2)
        a10 = (q1 - q2).dot(q0 - q2)
        a11 = (q1 - q2).dot(q1 - q2)
  
        A = np.matrix([[a00,a01],[a10,a11]])
        b = np.array([a00/2.0,a11/2.0])

        if A.determinant() == 0:

            candidate_indices.pop()
            candidate_center, candidate_indices = mcsc(points, candidate_indices)

        else:
            l = b*A.inverted()
  
            l_sum = l[0] + l[1]
            l2 = 1.0 - l_sum
        
            l_list = [l[0], l[1], l2]
  
            drop_index = None
            minimum = 0
  
            for index, number in enumerate(l_list):
                if number < minimum:
                    drop_index = index
                    minimum = number
     
            if drop_index != None:
                candidate_indices.pop(drop_index)
                candidate_center, candidate_indices = mcsc(points, candidate_indices)
            else:
                candidate_center = l[0]*q0 + l[1]*q1 + l2*q2  
 
    if candidate_count == 4:

        q0 = points[candidate_indices[0]]
        q1 = points[candidate_indices[1]]
        q2 = points[candidate_indices[2]]
        q3 = points[candidate_indices[3]]

        a00 = (q0 - q3).dot(q0 - q3)
        a01 = (q0 - q3).dot(q1 - q3)
        a02 = (q0 - q3).dot(q2 - q3)
        a10 = (q1 - q3).dot(q0 - q3)
        a11 = (q1 - q3).dot(q1 - q3)
        a12 = (q1 - q3).dot(q2 - q3)
        a20 = (q2 - q3).dot(q0 - q3)
        a21 = (q2 - q3).dot(q1 - q3)
        a22 = (q2 - q3).dot(q2 - q3)
  
        A = np.matrix([[a00, a01, a02],[a10, a11, a12],[a20, a21, a22]])
        b = np.array([a00/2.0, a11/2.0, a22/2.0])

        if A.determinant() == 0:

            candidate_indices.pop()
            candidate_center, candidate_indices = mcsc(points, candidate_indices)

        else:
            l = b*A.inverted()
  
            l_sum = l[0] + l[1] + l[2]
            l3 = 1.0 - l_sum
        
            l_list = [l[0], l[1], l[2], l3]
  
            drop_index = None
            minimum = 0
  
            for index, number in enumerate(l_list):
                if number < minimum:
                    drop_index = index
                    minimum = number
     
            if drop_index != None:
                candidate_indices.pop(drop_index)
                candidate_center, candidate_indices = mcsc(points, candidate_indices)

            else:
                candidate_center = l[0]*q0 + l[1]*q1 + l[2]*q2 + l3*q3  
        
    return candidate_center, candidate_indices

#reduce covering sphere along line p*(t - center)
#new candidates have p values between 0 exclusive and 1 exclusive
#candidate with smallest value is choosen for next covering sphere
#if minimum_covering sphere_candidate returns 4 points then minimum found
#if no new candidate points are found then minimum is found
#if new_center == center minimum is found

def find_next_candidate(points, center, candidate_indices):

    t = 0.0
    t, candidate_indices = mcsc(points, candidate_indices)

    if len(candidate_indices) == 4:
        return True, t, candidate_indices

    p = []
    for point in points:
        p.append(1)
  
    for index, point in enumerate(points):
        if index not in candidate_indices:
            d = (t - center).dot(points[candidate_indices[0]] - point)
            if d > 0:
                p[index] = (((points[candidate_indices[0]] + point)/2.0 - center)).dot( \
                            ((points[candidate_indices[0]] - point)/d))
    minimum = 1
    min_index = None
    for index, number in enumerate(p):
        if number < minimum and number > 0:
            min_index = index
            minimum = number
 
    if minimum == 1:
        return True, t, candidate_indices
 
    new_center = center + p[min_index]*(t - center)
 
    candidate_indices.insert(0, min_index) 
 
    if new_center.all() == center.all():
        return True, new_center, candidate_indices

    return False, new_center, candidate_indices

#given set of points defined
#return center and radius of minimum covering sphere
#define starting center as first point in points
#define p_1 as point farthest away from center
#first set of candidate points contains only p_1
#run until no new candidate points are found or until
#four candidate points are found that lie on minimum covering sphere

def minimum_covering_sphere(points):
        
    point_0 = points[0]
    new_center = point_0
    
    max_d = 0.0
    point_1_index = 0
    for index, point in enumerate(points):
        d = math.sqrt((point_0 - point).dot(point_0 - point))
        if d > max_d:
            point_1_index = index
            max_d = d
   
    candidate_indices = [point_1_index]
    finished = False

    while not finished:
        center = new_center
        finished, new_center, candidate_indices = \
        	find_next_candidate(points, center, candidate_indices)

    diff = points[candidate_indices[0]] - new_center
    radius = math.sqrt(diff.dot(diff))

    return new_center, radius


#calculate the minimum enclosing sphere volume of all elements (in um^3)
#uses the end points of filaments
#minimum covering sphere based on www.mel.nist.gov/msidlibrary/doc/hopp95.pdf
#returns a time and volume pair
def minEnclosingVolume(SnapshotList, snapshot=1):

	Snapshot = SnapshotList[i]

	points = []

	for filamentID, FS in Snapshot.filaments.iteritems():

		points.append(FS.coords[0])
		points.append(FS.coords[len(FS.coords) - 1])

	center, radius = minimum_covering_sphere(points)

	#return
	return (Snapshot.time, (4/3) * math.pi * pow(radius * 0.001,3))


#Calculate minimum enclosing sphere volumes from a number of trajectories (in um^3)
#uses the end points of filaments
#minimum covering sphere based on www.mel.nist.gov/msidlibrary/doc/hopp95.pdf
def minEnclosingVolumes(SnapshotLists, snapshot=1):

	minEnclosingVolumesList = []

	#get all data
	for SnapshotList in SnapshotLists:

		minEnclosingVolumes = minEnclosingVolume(SnapshotList, snapshot)
		minEnclosingVolumesList.append(minEnclosingVolumes)

	#now combine data
	finalTime = 0
	finalMinEnclosingVolume = 0

	for minEnclosingVolumes in minEnclosingVolumesList:

		finalTime += minEnclosingVolumes[0] 
		finalMinEnclosingVolume += minEnclosingVolumes[1]

	#average
	finalTime /= len(SnapshotLists)
	finalMinEnclosingVolume /= len(SnapshotLists)

	error = numpy.std(minEnclosingVolumesList)

	return (finalTime, finalMinEnclosingVolume, error)


#Calculate nematic order parameter of cylinders at each timestep
#from order parameter definition at http://cmt.dur.ac.uk/sjc/thesis_dlc/node65.html
#takes the largest eigenvalue of Q, a second rank ordering tensor over cylinder pairs
#returns a list of time and value pairs
def orderParameter(SnapshotList, snapshot=1, optype='apolar'):

	orderParams = []
	try:
		Snapshot = SnapshotList[snapshot]

	except IndexError:
		#print "Index error. Excluding this snapshot list."
		return False

	Q = np.zeros((3, 3))

	for alpha in [0,1,2]:

		for beta in [0,1,2]:

			for filamentID, FS in Snapshot.filaments.iteritems():

				#calculate direction
				mag = sqrt((FS.coords[len(FS.coords) - 1] - FS.coords[0]).dot(FS.coords[len(FS.coords) - 1] - FS.coords[0]))
				direction = (FS.coords[len(FS.coords) - 1] - FS.coords[0]) / mag

				#pairwise delta
				if alpha == beta: 
					delta = 1
				else:
					delta = 0

				#q value
				Q[alpha][beta] += ((3.0/2.0) * direction[alpha] * direction[beta] - (1.0/2.0) * delta)

			Q[alpha][beta] *= 1 / (float(len(Snapshot.filaments)))

	#take largest eigenvalue
	w, v= la.eig(Q)
	director = v[0]

	#apolar order parameter is largest eigenvalue
	if(optype == 'apolar'):

		return (Snapshot.time, max(w).real)

	#if we want polar order parameter
	elif(optype == 'polar'):

		sumCos = 0

		#now get cosines
		for filamentID, FS in Snapshot.filaments.iteritems():

			#calculate direction
			mag = sqrt((FS.coords[len(FS.coords) - 1] - FS.coords[0]).dot(FS.coords[len(FS.coords) - 1] - FS.coords[0]))
			direction = (FS.coords[len(FS.coords) - 1] - FS.coords[0]) / mag

			sumCos += director.dot(direction)

		sumCos *= 1 / (float(len(Snapshot.filaments)))

		return (Snapshot.time, abs(sumCos))


##Calculate nematic order parameter of cylinders from a number of trajectories
#from order parameter definition at http://cmt.dur.ac.uk/sjc/thesis_dlc/node65.html
#takes the largest eigenvalue of Q, a second rank ordering tensor over cylinder pairs
def orderParameters(SnapshotLists, snapshot=1, optype='apolar'):

	orderParamsList = []

	#get all data
	for SnapshotList in SnapshotLists:

		orderParam = orderParameter(SnapshotList, snapshot, optype)

		if(orderParam):
			orderParamsList.append(orderParam)

	finalTime = 0
	finalOrderParam = 0

	orderParamSet = []

	for orderParam in orderParamsList:

		finalTime += orderParam[0]
		finalOrderParam += orderParam[1]

		orderParamSet.append(orderParam[1])

	#average
	finalTime /= len(orderParamSet)
	finalOrderParam /= len(orderParamSet)

	error = numpy.std(orderParamSet)

	return (finalTime,finalOrderParam,error)

#Calculate the radius of gyration over all cylinders 
#returns a time and Rg pair for the chosen snapshot
def radiusOfGyration(SnapshotList, snapshot=1):

	Snapshot = SnapshotList[snapshot]

	#get center of mass of Snapshot
	com = np.array([0.0,0.0,0.0])
	beadCount = 0

		#go through filaments
	for filamentID, FS in Snapshot.filaments.iteritems():

		for coord in FS.coords:

			com += coord
			beadCount += 1

	#normalize
	com /= beadCount

	RgSquare = 0
	numCylinders = 0

	#Loop through cylinders
	for cylinderID, CS in Snapshot.cylinders.iteritems():

		numCylinders += 1
		cam = (CS.coords[1] + CS.coords[0]) / 2

		RgSquare += (cam - com).dot(cam - com)

	RgSquare /= numCylinders
	Rg = sqrt(RgSquare * 1e-6)

	return (Snapshot.time, Rg)


#Calculate the radius of gyration over all cylinders from a number of trajectories
def radiusOfGyrations(SnapshotLists, snapshot=1):

	RgList = []

	#get all data
	for SnapshotList in SnapshotLists:

		Rg = radiusOfGyration(SnapshotList, snapshot)
		RgList.append(Rg)

	#now combine data
	finalTime = 0
	finalRg = 0

	for Rg in RgList:

		finalTime += Rg[0]
		finalRg += Rg[1]

	#average
	finalTime /= len(RgList)
	finalRg /= len(RgList)

	Rgs = [x[1] for x in RgList]
	err = np.std(Rgs)

	return (finalTime,finalRg,err)

#Calculate contractile velocity towards the center of mass
#returns a list of time and velocity pairs
def contractileVelocity(SnapshotList):

	velocities = []

	for SnapshotIndex in xrange(0, len(SnapshotList) - 1):

		SnapshotI = SnapshotList[SnapshotIndex]
		SnapshotF = SnapshotList[SnapshotIndex + 1]

		deltaT = SnapshotF.time - SnapshotI.time

		#get center of mass of Snapshot i
		com = np.array([0.0,0.0,0.0])
		beadCount = 0

		#go through filaments
		for filamentID, FS in SnapshotI.filaments.iteritems():
			for coord in FS.coords:

				com += coord
				beadCount += 1

		#normalize
		com /= beadCount

		velocity = 0
		numCylinders = 0
		#Loop through cylinders
		for cylinderIDA, CSA in SnapshotF.cylinders.iteritems():

			#find in initial. if it doesnt exist, leave out
			try:
				CSB = SnapshotI.cylinders[cylinderIDA]
			except Exception, e:
				continue

			numCylinders += 1

			#calc velocity of midpoints of cylinders
			cam = (CSA.coords[1] + CSA.coords[0]) / 2
			cbm = (CSB.coords[1] + CSB.coords[0]) / 2

			vel = (cam - cbm)/ deltaT

			velocity += vel.dot((com - cbm)/((com - cbm).dot(com - cbm)))

		#normalize
		if numCylinders == 0: continue

		velocity /= numCylinders
		velocity *= 0.001 #units

		velocities.append((SnapshotF.time, velocity))

	#return velocities
	return velocities

#Calculate contractile velocity towards the center of mass from a number of trajectories
def contractileVelocities(SnapshotLists):

	velocitiesList = []

	#get all data
	for SnapshotList in SnapshotLists:

		velocities = contractileVelocity(SnapshotList)
		velocitiesList.append(velocities)

	#now combine data
	finalVelocities = []

	for i in xrange(0, len(velocities)):

		finalTime = 0
		finalVelocity = 0

		velocitySet = []

		for velocities in velocitiesList:

			finalTime += velocities[i][0]
			finalVelocity += velocities[i][1]

			velocitySet.append(velocities[i][1])

		#average
		finalTime /= len(SnapshotLists)
		finalVelocity /= len(SnapshotLists)

		error = numpy.std(velocitySet)

		finalVelocities.append((finalTime, finalVelocity, error))

	return finalVelocities

#Calculate average mean square displacement of filaments 
#returns a list of time and msd pairs
def meanSquareDisplacement(SnapshotList):

	msds = []

	#zero time point
	msds.append((0, 0))

	for SnapshotIndex in xrange(5, len(SnapshotList), 5):

		SnapshotI = SnapshotList[2]
		SnapshotF = SnapshotList[SnapshotIndex]

		msd = 0
		numFilaments = 0

		sds = []

		#Loop through cylinders
		for filamentIDA, FSA in SnapshotF.filaments.iteritems():

			#find in initial. if it doesnt exist, leave out
			try:
				FSB= SnapshotI.filaments[filamentIDA]
			except Exception, e:
				continue

			numFilaments += 1

			#calc com of fil A
			comA = np.array([0.0,0.0,0.0])
			beadCount = 0

			for coord in FSA.coords:

				comA += coord
				beadCount += 1

			#normalize
			comA /= beadCount

			#same for fil b
			comB = np.array([0.0,0.0,0.0])
			beadCount = 0

			for coord in FSB.coords:

				comB += coord
				beadCount += 1

			#normalize
			comB /= beadCount

			#calc square displacement (in um^2)
			sd = (comB - comA).dot(comB - comA) * 0.001 * 0.001

			sds.append(sd)

		#calc msd, error
		if numFilaments == 0: continue

		msd = np.mean(sds)
		msds.append((SnapshotF.time, msd))

	#return msds
	return msds

#Calculate average log mean square displacement of filaments over trajectories
#also returns log time for easy plotting
def meanSquareDisplacements(SnapshotLists):

	numTraj = 16
	numFilaments = 50

	msdsList = []

	#get all data
	for SnapshotList in SnapshotLists:

		msds = meanSquareDisplacement(SnapshotList)
		msdsList.append(msds)

	#now combine data
	finalMsds = []

	for i in xrange(0, len(msds)):

		finalTime = 0
		finalMsd = 0

		msdSet = []

		for msds in msdsList:

			finalTime += msds[i][0]
			finalMsd += msds[i][1]

			msdSet.append(msds[i][1])

		#average
		finalTime /= len(SnapshotLists)
		finalMsd /= len(SnapshotLists)

		#use standard error here
		err = numpy.std(msdSet) / sqrt(numTraj * numFilaments)

		finalMsds.append((finalTime, finalMsd, err))

	return finalMsds


#finds index of coordinate by linear mapping
def getIndex(coordinate, grid, compartment):

    index = 0
    i = 0

    for x in coordinate:
        if i == 0 :
            index += int(x / compartment[0])
        elif i == 1:
            index += int(x / compartment[1]) * grid[0];
        else:
            index += int(x / compartment[2]) * grid[0] * grid[1];
        i += 1

    return index


#Calculates density of filaments in each compartment for a Snapshot
#Uses number of cylinders in each compartment
#Param grid is the size of the grid (NX, NY, NZ)
#Param compartment is the size of compartment (nx, ny, nz)
#So, the entire grid is NX*nx * NY*ny * NZ*nz 
def density(Snapshot, grid, compartment):

	#create dict of compartment index and num cylinders by dimensions
	numCylinders = {}
	#create dict of compartment index and compartment
	compartmentMap = {}

	for i in xrange(0, grid[0]):

		for j in xrange(0, grid[1]):

			for k in xrange(0, grid[2]):

				coordinate = [i * compartment[0], \
							  j * compartment[1], \
							  k * compartment[2]]

				#add to map
				index = getIndex(coordinate, grid, compartment)
				compartmentMap[index] = coordinate

				#init num cylinders
				numCylinders[index] = 0

	#now, loop through cylinders, populate numCylinders
	for cylinderID, CS in Snapshot.cylinders.iteritems():

		#get coordinate, find index
		coordinate = (CS.coords[1] + CS.coords[0]) / 2
		index = getIndex(coordinate, grid, compartment)
	
		#add length of cylinder in um
		numCylinders[index] += 0.001 * sqrt((CS.coords[1] - CS.coords[0]). \
										 dot(CS.coords[1] - CS.coords[0]))

	#make final mapping of coordinates to num cylinders
	compartmentCoords = []
	densityValues = []

	for index, coordinate in compartmentMap.iteritems():

		compartmentCoords.append(compartmentMap[index])
		densityValues.append(numCylinders[index])

		##############################################
#
#
#
#			    PLOTTING ANALYSIS
#
#
#
###############################################

def calculateRgs(fileDirectory='', snapshot=1):

	FrameLists = readTrajectories('',fileDirectory=fileDirectory, runMin=0, runMax=16)
	return radiusOfGyrations(FrameLists, snapshot)

def calculateRgVsT(fileDirectory=''):

	FrameLists = readTrajectories('',fileDirectory=fileDirectory, runMin=0, runMax=15)

	Rgs1 = []

	for snap in xrange(2, len(FrameLists[0]), 2):
		Rgs1.append(radiusOfGyrations(FrameLists, snapshot=snap))

	del FrameLists

	return Rgs1

def calculateAllRgVsT():

	rgs = []

	rgs.append(calculateRgVsT(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR0.125/'))
	rgs.append(calculateRgVsT(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR0.25/'))
	rgs.append(calculateRgVsT(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR0.5/'))
	rgs.append(calculateRgVsT(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR2/'))
	rgs.append(calculateRgVsT(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR4/'))
	rgs.append(calculateRgVsT(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR8/'))

	return rgs

def calculateMsdVsT(fileDirectory=''):

	FrameLists = readTrajectories('',fileDirectory=fileDirectory, runMin=0, runMax=15)

	msds0 = meanSquareDisplacements(FrameLists)

	del FrameLists

	return msds0

def calculateAllMsdVsT():

	msds = []

	msds.append(calculateMsdVsT(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR0.125/'))
	msds.append(calculateMsdVsT(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR0.25/'))
	msds.append(calculateMsdVsT(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR0.5/'))
	msds.append(calculateMsdVsT(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/M0.02A0.1/'))
	msds.append(calculateMsdVsT(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR2/'))
	msds.append(calculateMsdVsT(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR4/'))
	msds.append(calculateMsdVsT(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR8/'))

	return msds

def calculateOrderParamVsT(fileDirectory='', numsnapshots=1, optype='apolar'):

	FrameLists = readTrajectories('',fileDirectory=fileDirectory, runMin=0, runMax=40)
	ops = []

	for i in xrange(0, len(FrameLists)):

		ops.append([])

		for j in xrange(0, numsnapshots):

			op = orderParameter(FrameLists[i], snapshot=j, optype=optype)

			if(op): 
				lastOp = op
				ops[i].append(op)
				di = 0

			#fill in data that did not run
			else:
				ops[i].append((lastOp[0] + 5.0*di, lastOp[1]))
				di += 1

	del FrameLists
	return ops

def calculateOrderParamPercentagesVsT(fileDirectory='', numsnapshots=1, optype='apolar'):

	FrameLists = readTrajectories('',fileDirectory=fileDirectory, runMin=0, runMax=40)
	ops = []

	for i in xrange(0, len(FrameLists)):

		ops.append([])

		for j in xrange(0, numsnapshots):

			op = orderParameter(FrameLists[i], snapshot=j, optype=optype)

			if(op): 
				lastOp = op
				ops[i].append(op)
				di = 0

			#fill in data that did not run
			else:
				ops[i].append((lastOp[0] + 5.0*di, lastOp[1]))
				di += 1

	del FrameLists

	#loop through runs, get percentages in range
	q1 = 0.0
	q2 = 0.33
	q3 = 0.66
	q4 = 1.0

	freqs = []

	for i in xrange(0, len(ops[0])):

		freqQ1 = 0
		freqQ2 = 0
		freqQ3 = 0

		time = 0

		for j in xrange(0, len(ops)):

			if(ops[j][i][1] >= q1 and ops[j][i][1] < q2):
				freqQ1 += 1

			if(ops[j][i][1] >= q2 and ops[j][i][1] < q3):
				freqQ2 += 1

			if(ops[j][i][1] >= q3 and ops[j][i][1] <= q4):
				freqQ3 += 1

			time += ops[j][i][0]

		freqs.append((time / len(ops), float(freqQ1) / len(ops), \
			float(freqQ2) / len(ops), float(freqQ3) / len(ops)))

	return freqs

def calculateOrderParamDistribution(fileDirectory='',numsnapshots=1, optype = 'apolar'):

	FrameLists = readTrajectories('',fileDirectory=fileDirectory, runMin=0, runMax=40)
	ops = []
	opsOneSnapshot = []

	for i in xrange(0, len(FrameLists)):

		ops.append([])

		for j in xrange(0, numsnapshots):

			op = orderParameter(FrameLists[i], snapshot=j, optype=optype)

			if(op): 
				lastOp = op
				ops[i].append(op)
				di = 0

			#fill in data that did not run
			else:
				ops[i].append((lastOp[0] + 5.0*di, lastOp[1]))
				di += 1

	del FrameLists
	for i in xrange(0, len(ops)):

		opsOneSnapshot.append(ops[i][snapshot][1])

	return opsOneSnapshot

def calculateOrderParamAverages(fileDirectory='', numsnapshots=1, optype='apolar'):

	FrameLists = readTrajectories('',fileDirectory=fileDirectory, runMin=0, runMax=20)
	ops = []

	for j in xrange(0, numsnapshots):

		ops.append(orderParameters(FrameLists, snapshot=j, optype=optype))

	del FrameLists

	return ops


def calculateOrderParamAverage(fileDirectory='', optype = 'polar', snapshot=1):

	FrameLists = readTrajectories('',fileDirectory=fileDirectory, runMin=0, runMax=15)

	op = orderParameters(FrameLists, snapshot=snapshot, optype=optype)

	del FrameLists

	return op

def calculateAllOrderParamAverage():

	ops = []

	ops.append(calculateOrderParamAverage(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR0.125/', optype='apolar', snapshot=200))
	ops.append(calculateOrderParamAverage(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR0.25/', optype='apolar', snapshot=200))
	ops.append(calculateOrderParamAverage(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR0.5/', optype='apolar', snapshot=200))
	ops.append(calculateOrderParamAverage(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/M0.02A0.1/', optype='apolar', snapshot=200))
	ops.append(calculateOrderParamAverage(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR2/', optype='apolar', snapshot=200))
	ops.append(calculateOrderParamAverage(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR4/', optype='apolar', snapshot=200))
	ops.append(calculateOrderParamAverage(fileDirectory='/Users/jameskomianos/Desktop/MEDYANData/medyandata/TR8/', optype='apolar', snapshot=200))

	return ops

def plotHistogram(histogramSnapshot, normalize=False):

	numBins = len(histogramSnapshot.hist)

	xfreqs = []
	xbins  = []

	#get bins and freqs
	index = 0
	for bin in sorted(histogramSnapshot.hist):

		binMin, binMax = bin

		xbins.append(binMin)
		xfreqs.append(histogramSnapshot.hist[bin])

		if(index == numBins):
			xbins.append(binMax)
		index += 1

	#normalize if needed
	if(normalize):
		maxfreq = max(xfreqs)
		xfreqs = [x/maxfreq for x in xfreqs]

	pos = np.arange(len(xbins))
	width = 0.1  # gives histogram aspect to the bar diagram

	ax = plt.axes()

	plt.bar(xbins, xfreqs, width, color='r')
	plt.show()

def plotRgHeatMap(saveFile=''):

	#read data
	matplotlib.rcParams['font.sans-serif']=["Arial"] 
	mpl.rcParams['xtick.labelsize'] = 8 
	mpl.rcParams['ytick.labelsize'] = 8 


	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	x_labels = [0.01,0.02,0.05,0.1,0.2,0.5]
	y_labels = [0.02, 0.01, 0.005]

	data = np.array([[1.460,1.392,1.199,1.144,1.105,1.120],
				     [1.600,1.518,1.387,1.332,1.305,1.268],
			         [1.713,1.652,1.512,1.437,1.449,1.410]])

	fig, ax = plt.subplots(figsize=(5,4))
	heatmap = ax.pcolor(data, cmap=plt.cm.YlGnBu_r)

	cbar = plt.colorbar(heatmap)

	cbar.ax.get_yaxis().set_ticks([])
	for j, lab in enumerate([1.0,1.2,1.4,1.6, 1.8]):
   		 cbar.ax.text(1.0, (1.5 * j) / 6.0, lab, ha='left', va='center', fontsize=9)

	cbar.ax.get_yaxis().labelpad = 30
	cbar.ax.set_ylabel(r'$\mathsf{R_g\/(\mu m)}$', fontsize=12)

	# put the major ticks at the middle of each cell
	ax.set_xticks(np.arange(data.shape[1])+0.5, minor=False)
	ax.set_yticks(np.arange(data.shape[0])+0.5, minor=False)

	# want a more natural, table-like display
	ax.invert_yaxis()

	ax.set_xticklabels(x_labels, minor=False)
	ax.set_yticklabels(y_labels, minor=False)
	plt.show()

	plt.ylabel(r'$\mathsf{R_{m:a}}$', fontsize=12)
	plt.xlabel(r'$\mathsf{R_{\alpha:a}}$', fontsize=12)

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)


def plotRgHeatMapExBind(saveFile=''):

	#read data
	x_labels = [0.01,0.02,0.05,0.1,0.2,0.5]
	y_labels = [0.02, 0.01, 0.005]

	data = np.array([[1.518,1.473,1.425,1.453,1.447,1.443],
				     [1.650,1.583,1.503,1.525,1.500,1.495],
			         [1.716,1.658,1.573,1.553,1.553,1.537]])

	fig, ax = plt.subplots()
	heatmap = ax.pcolor(data, cmap=plt.cm.YlGnBu_r)

	cbar = plt.colorbar(heatmap)

	cbar.ax.get_yaxis().set_ticks([])
	for j, lab in enumerate(['$1.4$','$1.5$','$1.6$','$1.7$','$1.8$']):
   		 cbar.ax.text(1.0, (1.5 * j) / 6.0, lab, ha='left', va='center')
	cbar.ax.get_yaxis().labelpad = 15
	cbar.ax.set_ylabel(r'$\mathrm{R_g\/(\mu m)}$', fontsize=16)


	# put the major ticks at the middle of each cell
	ax.set_xticks(np.arange(data.shape[1])+0.5, minor=False)
	ax.set_yticks(np.arange(data.shape[0])+0.5, minor=False)

	# want a more natural, table-like display
	ax.invert_yaxis()

	ax.set_xticklabels(x_labels, minor=False)
	ax.set_yticklabels(y_labels, minor=False)
	plt.show()

	plt.ylabel(r'$\mathrm{R_{m:a}}$', fontsize=20)
	plt.xlabel(r'$\mathrm{R_{\alpha:a}}$', fontsize=20)

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)

def plotOPHeatMap(saveFile=''):

	matplotlib.rcParams['font.sans-serif']=["Arial"] 
	mpl.rcParams['xtick.labelsize'] = 8 
	mpl.rcParams['ytick.labelsize'] = 8 

	#read data
	x_labels = [0.01,0.02,0.05,0.1,0.2,0.5]
	y_labels = [0.02, 0.01, 0.005]

	data = np.array([[0.144,0.275,0.527,0.618,0.691,0.637],
				     [0.153,0.145,0.313,0.389,0.477,0.451],
			         [0.133,0.139,0.188,0.241,0.276,0.321]])

	fig, ax = plt.subplots(figsize=(5,4))
	heatmap = ax.pcolor(data, cmap=plt.cm.YlGnBu_r)

	cbar = plt.colorbar(heatmap)

	cbar.ax.get_yaxis().set_ticks([])
	for j, lab in enumerate([0.1,0.2,0.3,0.4,0.5,0.6,0.7]):
   		 cbar.ax.text(1.2, (1.3 * j) / 8.0, lab, ha='left', va='center',fontsize=9)
	cbar.ax.get_yaxis().labelpad = 30
	cbar.ax.set_ylabel(r'$\mathsf{S}$', fontsize=8)

	# put the major ticks at the middle of each cell
	ax.set_xticks(np.arange(data.shape[1])+0.5, minor=False)
	ax.set_yticks(np.arange(data.shape[0])+0.5, minor=False)

	# want a more natural, table-like display
	ax.invert_yaxis()

	ax.set_xticklabels(x_labels, minor=False)
	ax.set_yticklabels(y_labels, minor=False)
	plt.show()

	plt.ylabel(r'$\mathsf{R_{m:a}}$', fontsize=8)
	plt.xlabel(r'$\mathsf{R_{\alpha:a}}$', fontsize=8)

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)


def plotRgRatioVsT(data, saveFile=''):

	numTraj = 16

	matplotlib.rcParams['font.sans-serif']=["Arial"] 
	mpl.rcParams['xtick.labelsize'] = 8 
	mpl.rcParams['ytick.labelsize'] = 8 

	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	fig, ax = plt.subplots(figsize=(3.3,3.3))

	#organize data
	data1 = data[0]
	data2 = data[1]
	data3 = data[2]
	data4 = data[3]
	data5 = data[4]
	data6 = data[5]

	time1 = [x[0] for x in data1]
	time2 = [x[0] for x in data2]
	time3 = [x[0] for x in data3]
	time4 = [x[0] for x in data4]
	time5 = [x[0] for x in data5]
	time6 = [x[0] for x in data6]

	Rg1_0 = data1[0][1]
	Rg2_0 = data2[0][1]
	Rg3_0 = data3[0][1]
	Rg4_0 = data4[0][1]
	Rg5_0 = data5[0][1]
	Rg6_0 = data6[0][1]

	err1_0 = data1[0][2]
	err2_0 = data2[0][2]
	err3_0 = data3[0][2]
	err4_0 = data4[0][2]
	err5_0 = data5[0][2]
	err6_0 = data6[0][2]

	Rg1 = [x[1]/Rg1_0 for x in data1]
	Rg2 = [x[1]/Rg2_0 for x in data2]
	Rg3 = [x[1]/Rg3_0 for x in data3]
	Rg4 = [x[1]/Rg4_0 for x in data4]
	Rg5 = [x[1]/Rg5_0 for x in data5]
	Rg6 = [x[1]/Rg6_0 for x in data6]

	#with propagation
	errplus_1 = [x[1]/Rg1_0 + (x[1]/Rg1_0)*sqrt(pow(x[2] / x[1], 2) + pow(err1_0 / Rg1_0, 2)) for x in data1]
	errminus_1 = [x[1]/Rg1_0 - (x[1]/Rg1_0)*sqrt(pow(x[2] / x[1], 2) + pow(err1_0 / Rg1_0, 2)) for x in data1]
	errplus_2 = [x[1]/Rg2_0 + (x[1]/Rg2_0)*sqrt(pow(x[2] / x[1], 2) + pow(err2_0 / Rg2_0, 2)) for x in data2]
	errminus_2 = [x[1]/Rg2_0 - (x[1]/Rg2_0)*sqrt(pow(x[2] / x[1], 2) + pow(err2_0 / Rg2_0, 2)) for x in data2]
	errplus_3 = [x[1]/Rg3_0 + (x[1]/Rg3_0)*sqrt(pow(x[2] / x[1], 2) + pow(err3_0 / Rg3_0, 2)) for x in data3]
	errminus_3 = [x[1]/Rg3_0 - (x[1]/Rg3_0)*sqrt(pow(x[2] / x[1], 2) + pow(err3_0 / Rg3_0, 2)) for x in data3]
	errplus_4 = [x[1]/Rg4_0 + (x[1]/Rg4_0)*sqrt(pow(x[2] / x[1], 2) + pow(err4_0 / Rg4_0, 2)) for x in data4]
	errminus_4 = [x[1]/Rg4_0 - (x[1]/Rg4_0)*sqrt(pow(x[2] / x[1], 2) + pow(err4_0 / Rg4_0, 2)) for x in data4]
	errplus_5 = [x[1]/Rg5_0 + (x[1]/Rg5_0)*sqrt(pow(x[2] / x[1], 2) + pow(err5_0 / Rg5_0, 2)) for x in data5]
	errminus_5 = [x[1]/Rg5_0 - (x[1]/Rg5_0)*sqrt(pow(x[2] / x[1], 2) + pow(err5_0 / Rg5_0, 2)) for x in data5]
	errplus_6 = [x[1]/Rg6_0 + (x[1]/Rg6_0)*sqrt(pow(x[2] / x[1], 2) + pow(err6_0 / Rg6_0, 2)) for x in data6]
	errminus_6 = [x[1]/Rg6_0 - (x[1]/Rg6_0)*sqrt(pow(x[2] / x[1], 2) + pow(err6_0 / Rg6_0, 2)) for x in data6]

	plt.plot(time1, Rg1, 'r', label=r'$\mathsf{\chi = 0.125}$', linewidth=2)
	plt.plot(time2, Rg2, 'g', label=r'$\mathsf{\chi = 0.25}$', linewidth=2)
	plt.plot(time3, Rg3, 'b', label=r'$\mathsf{\chi = 0.5}$', linewidth=2)
	plt.plot(time4, Rg4, 'y', label=r'$\mathsf{\chi = 2}$', linewidth=2)
	plt.plot(time5, Rg5, 'c', label=r'$\mathsf{\chi = 4}$', linewidth=2)
	plt.plot(time6, Rg6, 'm', label=r'$\mathsf{\chi = 8}$', linewidth=2)

	fill_between(time1, errplus_1, errminus_1, alpha=0.22, linewidth=2, color='r')
	fill_between(time2, errplus_2, errminus_2, alpha=0.22, linewidth=2, color='g')
	fill_between(time3, errplus_3, errminus_3, alpha=0.22, linewidth=2, color='b')
	fill_between(time4, errplus_4, errminus_4, alpha=0.22, linewidth=2, color='y')
	fill_between(time5, errplus_5, errminus_5, alpha=0.22, linewidth=2, color='c')
	fill_between(time6, errplus_6, errminus_6, alpha=0.22, linewidth=2, color='m')

	plt.xlabel(r'$\mathsf{Time\/(s)}$', fontsize=12)
	plt.ylabel(r'$\mathsf{R_{g,f} / R_{g,i}}$', fontsize=12)
	plt.xlim(0,2000)

	plt.legend(prop={'size':8})

	ax.yaxis.labelpad = 0
	ax.xaxis.labelpad = 0

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)


def plotMsdVsT(data, saveFile=''):

	matplotlib.rcParams['font.sans-serif']=["Arial"] 
	mpl.rcParams['xtick.labelsize'] = 8 
	mpl.rcParams['ytick.labelsize'] = 8 

	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	fig, ax = plt.subplots(figsize=(3.3,3.3))

	#organize data
	data1 = data[0]
	data2 = data[1]
	data3 = data[2]
	data4 = data[3]
	data5 = data[4]
	data6 = data[5]
	data7 = data[6]

	time1 = [x[0] for x in data1]
	time2 = [x[0] for x in data2]
	time3 = [x[0] for x in data3]
	time4 = [x[0] for x in data4]
	time5 = [x[0] for x in data5]
	time6 = [x[0] for x in data6]
	time7 = [x[0] for x in data7]

	msds1 = [x[1] for x in data1]
	msds2 = [x[1] for x in data2]
	msds3 = [x[1] for x in data3]
	msds4 = [x[1] for x in data4]
	msds5 = [x[1] for x in data5]
	msds6 = [x[1] for x in data6]
	msds7 = [x[1] for x in data7]

	plt.plot(time1, msds1,linewidth=4, color='r', label=r'$\mathsf{\chi = 0.125}$')
	plt.plot(time2, msds2,linewidth=4, color='g', label=r'$\mathsf{\chi = 0.25}$')
	plt.plot(time3, msds3,linewidth=4, color='b', label=r'$\mathsf{\chi = 0.5}$')
	plt.plot(time4, msds4,linewidth=4, color='k', label=r'$\mathsf{\chi = 1}$')
	plt.plot(time5, msds5,linewidth=4, color='y', label=r'$\mathsf{\chi = 2}$')
	plt.plot(time6, msds6,linewidth=4, color='c', label=r'$\mathsf{\chi = 4}$')
	plt.plot(time7, msds7,linewidth=4, color='m', label=r'$\mathsf{\chi = 8}$')

	plt.xlabel(r'$\mathsf{Time\/(s)}$', fontsize=12)
	plt.ylabel(r'$\mathsf{\langle \Delta x ^2 \rangle\/(\mu m^2)}$', fontsize=12)
	plt.xlim(0,2000)
	plt.legend(prop={'size':8})

	ax.yaxis.labelpad = 1
	ax.xaxis.labelpad = 0

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)


def plotDiffusionExponents(data, saveFile=''):

	matplotlib.rcParams['font.sans-serif']=["Arial"] 
	mpl.rcParams['xtick.labelsize'] = 8 
	mpl.rcParams['ytick.labelsize'] = 8 

	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	fig, ax = plt.subplots(figsize=(3.3,3.3))
	ax.get_xaxis().get_major_formatter().set_scientific(False)

	#only up to 1000s
	timeStart = 0
	timeCutoff = 1000

	#organize data
	data1 = data[0]
	data2 = data[1]
	data3 = data[2]
	data4 = data[3]
	data5 = data[4]
	data6 = data[5]
	data7 = data[6]

	time1 = [math.log(x[0], 10) for x in data1 if x[0] <= timeCutoff and x[0] != 0.0 and x[0] >= timeStart]
	msds1 = [math.log(x[1], 10) for x in data1 if x[0] <= timeCutoff and x[0] != 0.0 and x[0] >= timeStart]
	time2 = [math.log(x[0], 10) for x in data2 if x[0] <= timeCutoff and x[0] != 0.0 and x[0] >= timeStart]
	msds2 = [math.log(x[1], 10) for x in data2 if x[0] <= timeCutoff and x[0] != 0.0 and x[0] >= timeStart]
	time3 = [math.log(x[0], 10) for x in data3 if x[0] <= timeCutoff and x[0] != 0.0 and x[0] >= timeStart]
	msds3 = [math.log(x[1], 10) for x in data3 if x[0] <= timeCutoff and x[0] != 0.0 and x[0] >= timeStart]
	time4 = [math.log(x[0], 10) for x in data4 if x[0] <= timeCutoff and x[0] != 0.0 and x[0] >= timeStart]
	msds4 = [math.log(x[1], 10) for x in data4 if x[0] <= timeCutoff and x[0] != 0.0 and x[0] >= timeStart]
	time5 = [math.log(x[0], 10) for x in data5 if x[0] <= timeCutoff and x[0] != 0.0 and x[0] >= timeStart]
	msds5 = [math.log(x[1], 10) for x in data5 if x[0] <= timeCutoff and x[0] != 0.0 and x[0] >= timeStart]
	time6 = [math.log(x[0], 10) for x in data6 if x[0] <= timeCutoff and x[0] != 0.0 and x[0] >= timeStart]
	msds6 = [math.log(x[1], 10) for x in data6 if x[0] <= timeCutoff and x[0] != 0.0 and x[0] >= timeStart]
	time7 = [math.log(x[0], 10) for x in data7 if x[0] <= timeCutoff and x[0] != 0.0 and x[0] >= timeStart]
	msds7 = [math.log(x[1], 10) for x in data7 if x[0] <= timeCutoff and x[0] != 0.0 and x[0] >= timeStart]

	#linear regression for each
	treadmillingFactors = [0.125, 0.25, 0.5, 1, 2, 4, 8]
	exponents = []

	#get slopes
	slope, intercept, r_value, p_value, std_err = linregress(time1,msds1)
	exponents.append((slope, std_err))
	slope, intercept, r_value, p_value, std_err = linregress(time2,msds2)
	exponents.append((slope, std_err))
	slope, intercept, r_value, p_value, std_err = linregress(time3,msds3)
	exponents.append((slope, std_err))
	slope, intercept, r_value, p_value, std_err = linregress(time4,msds4)
	exponents.append((slope, std_err))
	slope, intercept, r_value, p_value, std_err = linregress(time5,msds5)
	exponents.append((slope, std_err))
	slope, intercept, r_value, p_value, std_err = linregress(time6,msds6)
	exponents.append((slope, std_err))
	slope, intercept, r_value, p_value, std_err = linregress(time7,msds7)
	exponents.append((slope, std_err))

	plt.errorbar(treadmillingFactors, [x[0] for x in exponents], yerr=[x[1] for x in exponents] ,fmt='o-', color='r')

	plt.xlabel(r'$\mathsf{\chi}$', fontsize=12)
	plt.ylabel(r'$\mathsf{\nu}}$', fontsize=12)
	#plt.xticks(np.arange(0.1, 10, 1.0))
	#plt.yticks(np.arange(0.0, 1.5, 0.25))

	ax.set_xscale('log')

	ax.yaxis.labelpad = 1
	ax.xaxis.labelpad = 0

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)

def plotOrderParamVsT(saveFile=''):

	data = calculateOrderParamVsT('')

	afont = {'fontname':'Arial'}

	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	#organize data
	for index in xrange(0, len(data)):

		data1 = data[index]

		time1 = [x[0] for x in data1]
		op1   = [x[1] for x in data1]

		plt.plot(time1, op1,linewidth=1, color='b')

	plt.xlabel('t', fontsize=24,  **afont)
	plt.ylabel('S', fontsize=24,  **afont)
	plt.xlim(0,2000)
	plt.legend()

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)


def plotOrderParamPercentagesVsT(saveFile=''):

	data = calculateOrderParamPercentagesVsT('/Users/jameskomianos/Desktop/Data/alignmentdata/Small2-MU2-NP', optype = 'polar')

	matplotlib.rcParams['font.sans-serif']=["Arial"] 

	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	#organize data
	colors = ['r', 'g', 'b']
	labels = [r'$\mathsf{0\/<\/\/S_1<\/1/3}$',r'$\mathsf{1/3\/<\/\/S_1<\/2/3}$', r'$\mathsf{2/3\/<\/\/S_1<\/1}$']

	for i in [0,1,2]:

		times = []
		percs = []

		for point in data:

			times.append(point[0])
			percs.append(point[1+i])

		plt.plot(times, percs, linewidth=1, color=colors[i],label=labels[i])

	plt.xlabel(r'$\mathsf{t}$', fontsize=24)
	plt.ylabel(r'$\mathsf{P(\Delta S_1)}$', fontsize=24)
	plt.xlim(0,2000)
	plt.ylim(0,1)
	plt.legend()

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)

def plotOrderParamDistribution(saveFile = ''):

	data = []

	#set data here
	data.append(calculateOrderParamDistribution('/Users/jameskomianos/Desktop/Data/alignmentdata/Small2-MU0.25-NP', snapshot=200, optype='apolar'))
	data.append(calculateOrderParamDistribution('/Users/jameskomianos/Desktop/Data/alignmentdata/Small2-MU0.25', snapshot=200, optype='apolar'))
	
	matplotlib.rcParams['font.sans-serif']=["Arial"] 

	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	#organize data
	colors = ['g', 'b']
	labels = [r'$\mathsf{xT}$',r'$\mathsf{T}$']

	numTrajs = 40.0

	for i in xrange(0, len(data)):

		data0 = data[i]

		normed_value = 1 / numTrajs

		hist, bins = np.histogram(data0, bins=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])

		widths = np.diff(bins) / len(data)
		hist = [x*normed_value for x in hist]

		plt.bar(bins[:-1] + (0.1/len(data)) * i, hist, widths, label=labels[i], color=colors[i], alpha=0.7)

	plt.xlabel(r'$\mathsf{S_2}$', fontsize=24)
	plt.ylabel(r'$\mathsf{P(S_2)}$', fontsize=24)

	plt.xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])

	plt.ylim(0,1)

	plt.legend()

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)	

def plotOrderParamBoxes(saveFile = ''):

	data = []

	#set data here
	data.append(calculateOrderParamDistribution('/Users/jameskomianos/Desktop/Data/alignmentdata/Small2-C2', snapshot=200, optype='apolar'))
	data.append(calculateOrderParamDistribution('/Users/jameskomianos/Desktop/Data/alignmentdata/Small2', snapshot=200, optype='apolar'))
	data.append(calculateOrderParamDistribution('/Users/jameskomianos/Desktop/Data/alignmentdata/Small2-C22', snapshot=200, optype='apolar'))

	matplotlib.rcParams['font.sans-serif']=["Arial"] 

	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	#organize data
	plt.boxplot([data[0],data[1],data[2]], sym='')

	labels = [2,12,22]
	plt.xticks([1,2,3], labels)

	plt.ylabel(r'$\mathsf{S_2}$', fontsize=24)
	plt.xlabel(r'$\mathsf{F_0 (pN)}$', fontsize=24)

	plt.legend()

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)

def plotOrderParamAverages(saveFile = ''):

	data = []

	#set data here
	data.append(calculateOrderParamAverages('/Users/jameskomianos/Desktop/Data/alignmentdata/Small4', optype='polar'))
	data.append(calculateOrderParamAverages('/Users/jameskomianos/Desktop/Data/alignmentdata/Small8', optype='polar'))
	data.append(calculateOrderParamAverages('/Users/jameskomianos/Desktop/Data/alignmentdata/Small16', optype='polar'))

	print data

	matplotlib.rcParams['font.sans-serif']=["Arial"] 

	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	colors = ['g', 'b', 'r']
	labels = [r'$\mathsf{N=4}$',r'$\mathsf{N=8}$', r'$\mathsf{N=16}$']

	#organize data
	for index in xrange(0, len(data)):

		data1 = data[index]

		time1 = [x[0] for x in data1]
		op1   = [x[1] for x in data1]
		err1  = [x[2] for x in data1]

		plt.errorbar(time1, op1, yerr=err1, linewidth=1, color=colors[index], label = labels[index])

	plt.xlabel('t', fontsize=24)
	plt.ylabel(r'$\mathsf{S_1}$', fontsize=24)
	plt.xlim(0,2000)
	plt.legend()

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)

def plotOPHeatMap2(saveFile=''):

	matplotlib.rcParams['font.sans-serif']=["Arial"] 

	#read data
	x_labels = [10,40,100]
	y_labels = [1.0, 0.25, 0.1]

	data = np.array([[0.354,0.124,0.215],
				     [0.263,0.218,0.207],
			         [0.223,0.193,0.161]])

	fig, ax = plt.subplots()
	heatmap = ax.pcolor(data, cmap=plt.cm.RdBu_r, vmin=0, vmax=1)

	cbar = plt.colorbar(heatmap)

	cbar.set_clim(0, 1.0)

	cbar.ax.get_yaxis().set_ticks([])
	for j, lab in enumerate([r'$\mathsf{0}$',r'$\mathsf{0.2}$',r'$\mathsf{0.4}$',r'$\mathsf{0.6}$',r'$\mathsf{0.8}$', r'$\mathsf{1.0}$']):
   		 cbar.ax.text(1.0, (1.4 * j) / 7.0, lab, ha='left', va='center')
	cbar.ax.get_yaxis().labelpad = 20
	cbar.ax.set_ylabel(r'$\mathsf{S_1}$', fontsize=24)

	# put the major ticks at the middle of each cell
	ax.set_xticks(np.arange(data.shape[1])+0.5, minor=False)
	ax.set_yticks(np.arange(data.shape[0])+0.5, minor=False)

	# want a more natural, table-like display
	ax.invert_yaxis()

	ax.set_xticklabels(x_labels, minor=False)
	ax.set_yticklabels(y_labels, minor=False)
	plt.show()

	plt.ylabel(r'$\mathsf{k_u (1 / s)}$', fontsize=20)
	plt.xlabel(r'$\mathsf{V_m (nm / s)}$', fontsize=20)

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)


def plotOrderParameterVsTreadmill(data, saveFile=''):

	matplotlib.rcParams['font.sans-serif']=["Arial"] 
	mpl.rcParams['xtick.labelsize'] = 8 
	mpl.rcParams['ytick.labelsize'] = 8 

	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	fig, ax = plt.subplots(figsize=(3.3,3.3))
	ax.get_xaxis().get_major_formatter().set_scientific(False)

	#organize data
	ops = [x[1] for x in data]
	err  = [x[2] for x in data]
	treadmillingFactors = [0.125, 0.25, 0.5, 1, 2, 4, 8]

	plt.errorbar(treadmillingFactors, ops, yerr=err, fmt='o-', color='r')

	plt.xlabel(r'$\mathsf{\chi}$', fontsize=12)
	plt.ylabel(r'$\mathsf{\nu}$', fontsize=12)

	ax.set_xscale('log')

	ax.yaxis.labelpad = 5
	ax.xaxis.labelpad = 0

	#plt.yticks(np.arange(0.2, 0.75, 0.1))

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)
