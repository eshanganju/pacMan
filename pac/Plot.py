""" Potting module
"""


import os
import sys
import numpy
import random
import math

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from pac import Measure

# This is to plot all the text when running the functions
VERBOSE = True
TESTING = True

def plotPSD(gsd0, gsd1, gsd2=np.zeros((1,2)), gsd3=np.zeros((1,2)), gsd4=np.zeros((1,2)), gsd5=np.zeros((1,2)), gsd6=np.zeros((1,2)), xmax = 10, xmin = 0.001):
	"""
	plots semilogx grain size distribution
	add additional  gradations if passed
	add a legend with the gradations, ask user for legend labels
	"""
	plt.figure()
	plt.semilogx( gsd0[ :,0 ] , gsd0[ : , 1 ] , label='0' )
	plt.semilogx( gsd1[ :,0 ] , gsd1[ : , 1 ] , label='1' )
	if gsd2.sum() != 0: plt.semilogx( gsd2[ :,0 ] , gsd2[ : , 1 ] , label='2' )
	if gsd3.sum() != 0: plt.semilogx( gsd3[ :,0 ] , gsd3[ : , 1 ] , label='3' )
	if gsd4.sum() != 0: plt.semilogx( gsd4[ :,0 ] , gsd4[ : , 1 ] , label='4' )
	if gsd5.sum() != 0: plt.semilogx( gsd5[ :,0 ] , gsd5[ : , 1 ] , label='5' )
	if gsd6.sum() != 0: plt.semilogx( gsd6[ :,0 ] , gsd6[ : , 1 ] , label='6' )
	plt.xlim( 0.01 , 10 )
	plt.ylim( 0 , 100 )
	plt.legend()
	plt.show()


def rosePlot():
	"""
	"""


def blobPlot():
	"""
	"""


def plotOrientationsSPAM(orientations_zyx,
							projection="lambert",
							plot="both",
							binValueMin=None,
							binValueMax=None,
							binNormalisation = False,
							numberOfRings = 9,
							pointMarkerSize = 8,
							cmap = matplotlib.pyplot.cm.RdBu_r,
							title = "",
							subtitle = {"points":"","bins":""},
							saveFigPath = None,
							sampleName='',
							figXSize = 6.1,
							figYSize = 4.8,
							figFontSize = 15,
							labelName = '' ):
	"""
	2021-04-21: This function is taken from the SPAM library. Modifications have been made
	the asthetics of the figure to have the figures in publication quality. 

	Majority of the code is kept the same.

	Main function for plotting 3D orientations.
	This function plots orientations (described by unit-direction vectors) from a sphere onto a plane.

	One useful trick for evaluating these orientations is to project them with a "Lambert equal area projection", which means that an isotropic distribution of angles is projected as equally filling the projected space.

	Parameters
	----------
		orientations : Nx3 numpy array of floats
			Z, Y and X components of direction vectors.
			Non-unit vectors are normalised.

		projection : string, optional
			Selects different projection modes:
				**lambert** : Equal-area projection, default and highly reccommended. See https://en.wikipedia.org/wiki/Lambert_azimuthal_equal-area_projection

				**equidistant** : equidistant projection

		plot : string, optional
			Selects which plots to show:
				**points** : shows projected points individually
				**bins** : shows binned orientations with counts inside each bin as colour
				**both** : shows both representations side-by-side, default

		title : string, optional
			Plot main title. Default = ""

		subtitle : dictionary, optional
			Sub-plot titles:
				**points** : Title for points plot. Default = ""
				**bins** : Title for bins plot. Default = ""

		binValueMin : int, optional
			Minimum colour-bar limits for bin view.
			Default = None (`i.e.`, auto-set)

		binValueMax : int, optional
			Maxmum colour-bar limits for bin view.
			Default = None (`i.e.`, auto-set)

		binNormalisation : bool, optional
			In binning mode, should bin counts be normalised by mean counts on all bins
			or absoulte counts?

		cmap : matplotlib colour map, optional
			Colourmap for number of counts in each bin in the bin view.
			Default = ``matplotlib.pyplot.cm.RdBu_r``

		numberOfRings : int, optional
			Number of rings (`i.e.`, radial bins) for the bin view.
			The other bins are set automatically to have uniform sized bins using an algorithm from Jacquet and Tabot.
			Default = 9 (quite small bins)

		pointMarkerSize : int, optional
			Size of points in point view.
			Default = 8 (quite big points)

		saveFigPath : string, optional
			Path to save figure to -- stops the graphical plotting.
			Default = None

	Returns
	-------
		None -- A matplotlib graph is created and show()n

	Note
	----
		Authors: Edward Andò, Hugues Talbot, Clara Jacquet and Max Wiebicke
	"""
	import matplotlib.pyplot as plt
	# ========================================================================
	# ==== Reading in data, and formatting to x,y,z sphere                 ===
	# ========================================================================
	numberOfPoints = orientations_zyx.shape[0]

	# ========================================================================
	# ==== Check that all the vectors are unit vectors                     ===
	# ========================================================================
	if VERBOSE: print( "\t-> Normalising all vectors in x-y-z representation..." ),

	# from http://stackoverflow.com/questions/2850743/numpy-how-to-quickly-normalize-many-vectors
	norms = numpy.apply_along_axis( numpy.linalg.norm, 1, orientations_zyx )
	orientations_zyx = orientations_zyx / norms.reshape( -1, 1 )

	if VERBOSE: print( "done." )

	# ========================================================================
	# ==== At this point we should have clean x,y,z data in memory         ===
	# ========================================================================
	if VERBOSE: print( "\t-> We have %i orientations in memory."%( numberOfPoints ) )

	projection_xy       = numpy.zeros( (numberOfPoints, 2) )

	projection_theta_r  = numpy.zeros( (numberOfPoints, 2) )

	# ========================================================================
	# ==== Projecting from x,y,z sphere to the desired projection          ===
	# ========================================================================
	for vectorN in range( numberOfPoints ):
		z,y,x = orientations_zyx[ vectorN ]

		if z < 0: z = -z; y = -y; x = -x

		projection_xy[ vectorN ], projection_theta_r[ vectorN ] = _projectOrientations( projection, "xyz", [x,y,z] )

	if projection == "lambert": radiusMax = numpy.sqrt(2)
	elif projection == "stereo": radiusMax = 1.0
	elif projection == "equidistant": radiusMax = 1.0

	if VERBOSE: print( "\t-> Biggest projected radius (r,t) = {}".format( numpy.abs( projection_theta_r[:,1] ).max() ) )


	if plot == "points" or plot == "both":
		fig = plt.figure( figsize = ( figXSize, figYSize ) )
		fig.suptitle( title )

		if plot == "both":
			ax  = fig.add_subplot( 121, polar=True )

		else:
			ax  = fig.add_subplot( 111, polar=True )

		ax.set_title( subtitle['points']+"\n" )
		plt.axis( ( 0, math.pi*2, 0, radiusMax ) )
		plt.xticks(fontsize=figFontSize)

		radiusGridAngles = numpy.arange( 15, 91, 15 )
		radiusGridValues = []

		for angle in radiusGridAngles:
			radiusGridValues.append( _projectOrientations( projection, "spherical", [ 1, angle*math.pi/180.0, 0 ] )[1][1] )

		ax.set_rgrids( radiusGridValues, labels=[ "%02i$^\circ$"%(x) for x in numpy.arange(  15,91,15) ], angle=None, fmt=None, fontsize=12)
		ax.plot( projection_theta_r[:,0], projection_theta_r[:,1] , '.', markersize=pointMarkerSize )

		if plot == "points":
			plt.show()

	if plot == "bins" or plot == "both":
		import matplotlib.patches
		import matplotlib.colorbar
		import matplotlib.collections

		if plot == "both":
			ax  = fig.add_subplot( 122, polar=True )

		if plot == "bins":
			fig = plt.figure()
			ax  = fig.add_subplot( 111, polar=True)

		if VERBOSE: print( "\t-> Starting Data binning..." )

		if VERBOSE: print( "\t-> Number of Rings (radial bins) = ", numberOfRings )

		numberOfAngularBinsPerRing = numpy.arange( 1, numberOfRings+1, 1 )
		numberOfAngularBinsPerRing = 4 * ( 2 * numberOfAngularBinsPerRing - 1 )

		if VERBOSE: print( "\t-> Number of angular bins per ring = ", numberOfAngularBinsPerRing )

		binCounts = numpy.zeros( ( numberOfRings, numberOfAngularBinsPerRing[-1] ) )

		# ========================================================================
		# ==== Start counting the vectors into bins                            ===
		# ========================================================================
		for vectorN in range( numberOfPoints ):
			angle, radius = projection_theta_r[ vectorN, : ]

			if angle < 0:             angle += 2*math.pi
			if angle > 2 * math.pi:   angle -= 2*math.pi

			ringNumber = int(numpy.floor( radius / ( radiusMax / float(numberOfRings) ) ) )

			if ringNumber > numberOfRings - 1:
				if VERBOSE: print( "\t-> Point with projected radius = %f is a problem (radiusMax = %f), putting in furthest  bin"%( radius, radiusMax ) )
				ringNumber = numberOfRings - 1

			angularBin = int( numpy.floor( ( angle ) / ( 2 * math.pi / float( numberOfAngularBinsPerRing[ ringNumber ] ) ) ) ) + 1

			if angularBin > numberOfAngularBinsPerRing[ringNumber] - 1:
				if VERBOSE: print( "\t-> Point with projected angle = %f does not belong to the last bin, putting in first bin"%( angle ) )
				angularBin = 0

			binCounts[ ringNumber, angularBin ] += 1

		# ========================================================================
		# === Plotting binned data                                             ===
		# ========================================================================

		plottingRadii = numpy.linspace( radiusMax/float(numberOfRings), radiusMax, numberOfRings )
		bars = []

		'''
		add two fake, small circles at the beginning so that they are overwritten
		they will be coloured with the min and max colour
		theta   radius    width
		'''

		bars.append( [   0,   radiusMax,    2*math.pi ] )
		bars.append( [   0,   radiusMax,    2*math.pi ] )

		flatBinCounts = numpy.zeros( numpy.sum( numberOfAngularBinsPerRing ) + 2 )

		binNumber = 2

		if binNormalisation:
			avg_binCount = float(numberOfPoints)/numpy.sum( numberOfAngularBinsPerRing )
			if VERBOSE: print( "\t-> Average binCount = ", avg_binCount )

		for ringNumber in range( numberOfRings )[::-1]:
			deltaTheta    = 360 / float( numberOfAngularBinsPerRing[ringNumber] )
			deltaThetaRad = 2 * math.pi / float( numberOfAngularBinsPerRing[ringNumber] )

			# --- Angular bins                                                 ---
			for angularBin in range( numberOfAngularBinsPerRing[ringNumber] ):

				bars.append( [ angularBin*deltaThetaRad - deltaThetaRad/2.0, plottingRadii[ ringNumber ], deltaThetaRad ] )
				#bars.append( ax.bar( angularBin*deltaThetaRad - deltaThetaRad/2.0, plottingRadii[ ringNumber ], deltaThetaRad, bottom=0.0 ) )

				# Add the number of vectors counted for this bin
				if binNormalisation:
					flatBinCounts[ binNumber ] = binCounts[ ringNumber, angularBin ]/avg_binCount
				else:
					flatBinCounts[ binNumber ] = binCounts[ ringNumber, angularBin ]

				# Add one to bin number
				binNumber += 1

		del binNumber

		# figure out auto values if they're requested.
		if binValueMin is None: binValueMin = flatBinCounts[2::].min()
		if binValueMax is None: binValueMax = flatBinCounts[2::].max()

		# Add two flat values for the initial wedges.
		flatBinCounts[0] = binValueMin
		flatBinCounts[1] = binValueMax

		barsPlot = ax.bar( numpy.array( bars )[:,0], numpy.array( bars )[:,1], width=numpy.array( bars )[:,2], bottom=0.0)

		for binCount,bar in zip(  flatBinCounts, barsPlot ):
			bar.set_facecolor( cmap(  ( binCount - binValueMin) / float( binValueMax - binValueMin ) ) )

		plt.axis( [ 0, numpy.deg2rad(360), 0, radiusMax ] )
		plt.xticks(fontsize=figFontSize)

		norm = matplotlib.colors.Normalize( vmin=binValueMin, vmax=binValueMax )

		ax3 = fig.add_axes([0.9, 0.1, 0.03, 0.8])
		cb1 = matplotlib.colorbar.ColorbarBase( ax3, cmap=cmap, norm=norm, ticks=[0,2,4,6,8,10,12,14,16,18,20])
		cb1.set_label(label=labelName, size=figFontSize)
		cb1.ax.tick_params(labelsize=figFontSize)

		radiusGridAngles = numpy.arange( 15, 91, 15 )
		radiusGridValues = []

		for angle in radiusGridAngles:
			radiusGridValues.append( _projectOrientations( projection, "spherical", [ 1, angle*math.pi/180.0, 0 ] )[1][1] )

		ax.set_rgrids( radiusGridValues, labels=[ "%02i$^\circ$"%(x) for x in numpy.arange(  15,91,15) ], angle=None, fmt=None, fontsize=12)

		fig.subplots_adjust(left=0.05,right=0.85)

		plt.rcParams['font.size'] = figFontSize

		if saveFigPath is not None:
			plt.draw()
			plt.savefig( (saveFigPath+sampleName), transparent=False, dpi=300 )
			plt.pause(1)
			plt.close()

		else:
			plt.draw()
			plt.pause(1)
			plt.close()


def _projectOrientations(projection, coordSystem, vector):
	"""1. type of projection
	2. coordinate system (XYZ, or spherical)
	3. input 3-unit vector, or array:
		for x-y-z that order is maintained
		for spherical: r, theta (inclination), phi (azimuth)
	"""
	# Returns: projection_xy, projection_theta_r

	projection_xy_local = numpy.zeros( (2) )
	projection_theta_r_local = numpy.zeros( (2) )

	if coordSystem == "spherical":
		# unpack vector 
		r, theta, phi = vector

		x = r * math.sin( theta ) * math.cos( phi )
		y = r * math.sin( theta ) * math.sin( phi )
		z = r * math.cos( theta )

	else:
		# unpack vector 
		x, y, z = vector
		# we're in cartesian coordinates, (x-y-z mode) Calculate spherical coordinates
		# passing to 3d spherical coordinates too...
		# From: https://en.wikipedia.org/wiki/Spherical_coordinate_system
		#  Several different conventions exist for representing the three coordinates, and for the order in which they should be written.
		#  The use of (r, θ, φ) to denote radial distance, inclination (or elevation), and azimuth, respectively, is common practice in physics, 
		#   and is specified by ISO standard 80000-2 :2009, and earlier in ISO 31-11 (1992).
		r = numpy.sqrt( x**2 + y**2 + z**2 )
		theta = math.acos( z / r )   # inclination
		phi = math.atan2( y, x )   # azimuth

	if projection == "lambert":			# dividing by sqrt(2) so that we're projecting onto a unit circle
		projection_xy_local[0] = x*( math.sqrt(2/(1+z)) )
		projection_xy_local[1] = y*( math.sqrt(2/(1+z)) )

		# sperhical coordinates -- CAREFUL as per this wikipedia page: https://en.wikipedia.org/wiki/Lambert_azimuthal_equal-area_projection
		# the symbols for inclination and azimuth ARE INVERTED WITH RESPEST TO THE SPHERICAL COORDS!!!
		projection_theta_r_local[0] = phi
		# HACK: doing math.pi - angle in order for the +z to be projected to 0,0
		projection_theta_r_local[1] = 2 * math.cos( ( math.pi - theta ) / 2 )

		# cylindrical coordinates
		#projection_theta_r_local[0] = phi
		#projection_theta_r_local[1] = math.sqrt( 2.0 * ( 1 + z ) )

	if projection == "stereo":
		projection_xy_local[0] = x / ( 1 - z )
		projection_xy_local[1] = y / ( 1 - z )

		# https://en.wikipedia.org/wiki/Stereographic_projection uses a different standard from the page on spherical coord Spherical_coordinate_system
		projection_theta_r_local[0] = phi
		# HACK: doing math.pi - angle in order for the +z to be projected to 0,0
		# HACK: doing math.pi - angle in order for the +z to be projected to 0,0
		projection_theta_r_local[1] = numpy.sin( math.pi - theta ) / ( 1 - numpy.cos( math.pi - theta ) )

	if projection == "equidistant":
		# https://en.wikipedia.org/wiki/Azimuthal_equidistant_projection
		# TODO: To be checked, but this looks like it should -- a straight down projection.
		projection_xy_local[0] = math.sin( phi )
		projection_xy_local[1] = math.cos( phi )

		projection_theta_r_local[0] = phi
		projection_theta_r_local[1] = numpy.cos( theta - math.pi/2 )

	return projection_xy_local, projection_theta_r_local


def plotLabelPropertyMap( lableMap, propertyToPlot='volume', sliceNumber=0,
							sampleName='',saveImg=True,outputDir=''):
	"""Plot a slice of the tomography and color the particles with the chosen property
	"""
	print('pffft')