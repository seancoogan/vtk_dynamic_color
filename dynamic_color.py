#!/usr/bin/env python
 
from __future__ import print_function
import math
import vtk
import sys
import random
from random import randint
from color import Color

# The constant values for the colorblind options
ALL = 1
DEUTERANOPIA = 2
PROTONOPIA = 3
TRITANOPIA = 4
isColorblind = 0 # start by assuming not colorblind
total_bands = 8

# 
# True: use Grey as the midpoint
# False: use the midpoint calculation algorithm
#
GREY_VALUE = True

# The Display Colors
dispColor1 = Color()
dispColor2 = Color()

'''
	Divide a range into bands
	:param: dR - [min, max] the range that is to be covered by the bands.
	:param: numberOfBands - the number of bands, a positive integer.
	:param: nearestInteger - if True then [floor(min), ceil(max)] is used.
	:return: A List consisting of [min, midpoint, max] for each band.
'''
def MakeBands(dR, numberOfBands, nearestInteger):
	bands = list()
	if (dR[1] < dR[0]) or (numberOfBands <= 0):
		return bands
	x = list(dR)
	if nearestInteger:
		x[0] = math.floor(x[0])
		x[1] = math.ceil(x[1])
	dx = (x[1] - x[0])/float(numberOfBands)
	b = [x[0], x[0] + dx / 2.0, x[0] + dx]
	i = 0
	while i < numberOfBands:
		bands.append(b)
		b = [b[0] + dx, b[1] + dx, b[2] + dx]
		i += 1
	return bands

'''
	Divide a range into integral bands
	:param: dR - [min, max] the range that is to be covered by the bands.
	:return: A List consisting of [min, midpoint, max] for each band.
'''
def MakeIntegralBands(dR):
	bands = list()
	if (dR[1] < dR[0]):
		return bands
	x = list(dR)
	x[0] = math.floor(x[0])
	x[1] = math.ceil(x[1])
	numberOfBands = int(abs(x[1]) + abs(x[0]))
	return MakeBands(x,numberOfBands, False)
	
'''
	Generate elevations over the surface.
	:param: src - the vtkPolyData source.
	:return: - vtkPolyData source with elevations.
'''
def MakeElevations(src):

	bounds = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
	src.GetBounds(bounds)
	elevFilter = vtk.vtkElevationFilter()
	elevFilter.SetInputData(src)
	elevFilter.SetLowPoint(0, bounds[2], 0)
	elevFilter.SetHighPoint(0, bounds[3], 0)
	elevFilter.SetScalarRange(bounds[2], bounds[3])
	elevFilter.Update()
	return elevFilter.GetPolyDataOutput()
	
'''
	Make a parametric surface as the source.
	:return: vtkPolyData with normal and scalar data.
'''
def MakeParametricSource():

	fn = vtk.vtkParametricRandomHills()
	fn.AllowRandomGenerationOn()
	fn.SetRandomSeed(1)
	fn.SetNumberOfHills(30)
	if fn.GetClassName() == 'vtkParametricRandomHills':
		# Make the normals face out of the surface.
		fn.ClockwiseOrderingOff()

	source = vtk.vtkParametricFunctionSource()
	source.SetParametricFunction(fn)
	source.SetUResolution(50)
	source.SetVResolution(50)
	source.SetScalarModeToZ()
	source.Update()
	# Name the arrays (not needed in VTK 6.2+ for vtkParametricFunctionSource)
	source.GetOutput().GetPointData().GetNormals().SetName('Normals')
	source.GetOutput().GetPointData().GetScalars().SetName('Scalars')
	return source.GetOutput()

'''
	Make a lookup table using vtkColorSeries.
	:return: An indexed lookup table.
'''
def MakeLUT(num_distinct):
	
	if num_distinct < 3:
		print("The number of bands must be greater than or equal to 3.")
		sys.exit(1)
	else:
		ctf = vtk.vtkColorTransferFunction()
		ctf.SetColorSpaceToDiverging()
		
		ctf.AddRGBPoint(0.0, dispColor1.r, dispColor1.g, dispColor1.b)
		
		if GREY_VALUE:
			ctf.AddRGBPoint(0.5, 0.5, 0.5, 0.5) # grey value 
		else:
			midpoint = dispColor1.calculate_midpoint(dispColor2)
			ctf.AddRGBPoint(0.5, midpoint[0], midpoint[1], midpoint[2])
			
		ctf.AddRGBPoint(1.0, dispColor2.r, dispColor2.g, dispColor2.b)
 
		lut = vtk.vtkLookupTable()
		lut.SetNumberOfTableValues(num_distinct)
		lut.Build()
    

		for i in range(0, num_distinct):
			rgb = list(ctf.GetColor(float(i)/(num_distinct+1)))+[1]
			
			# set the value into the table
			lut.SetTableValue(i,rgb)
			
	return lut

'''
	Create a lookup table with the colors reversed.
	:param: lut - An indexed lookup table.
	:return: The reversed indexed lookup table.
'''
def ReverseLUT(lut):

	lutr = vtk.vtkLookupTable()
	lutr.DeepCopy(lut)
	t = lut.GetNumberOfTableValues() - 1
	revList = reversed(list(range(t + 1)))
	for i in revList:
		rgba = [0,0,0]
		v = float(i)
		lut.GetColor(v,rgba)
		rgba.append(lut.GetOpacity(v))
		lutr.SetTableValue(t - i,rgba)
	t = lut.GetNumberOfAnnotatedValues() - 1
	for i in revList:
		lutr.SetAnnotation(t - i, lut.GetAnnotation(i))
	return lutr

'''
    Count the number of scalars in each band.
    :param: bands - the bands.
    :param: src - the vtkPolyData source.
    :return: The frequencies of the scalars in each band.
'''
def Frequencies(bands, src):

    freq = dict()
    for i in range(len(bands)):
        freq[i] = 0;
    tuples = src.GetPointData().GetScalars().GetNumberOfTuples()
    for i in range(tuples):
        x = src.GetPointData().GetScalars().GetTuple1(i)
        for j in range(len(bands)):
            if x <= bands[j][2]:
                freq[j] = freq[j] + 1
                break
    return freq

'''
	Glyph the normals on the surface.
	You may need to adjust the parameters for maskPts, 
	arrow and glyph for a nice appearance.
 
	:param: src - the surface to glyph.
	:param: reverseNormals - if True the normals on the surface 
		are reversed.
	:return: The glyph object.
'''
def MakeGlyphs(src, reverseNormals):

	# Sometimes the contouring algorithm can create a volume whose gradient
	# vector and ordering of polygon (using the right hand rule) are
	# inconsistent. vtkReverseSense cures this problem.
	reverse = vtk.vtkReverseSense()
 
	# Choose a random subset of points.
	maskPts = vtk.vtkMaskPoints()
	maskPts.SetOnRatio(5)
	maskPts.RandomModeOn()
	if reverseNormals:
		reverse.SetInputData(src)
		reverse.ReverseCellsOn()
		reverse.ReverseNormalsOn()
		maskPts.SetInputConnection(reverse.GetOutputPort())
	else:
		maskPts.SetInputData(src)
 
	# Source for the glyph filter
	arrow = vtk.vtkArrowSource()
	arrow.SetTipResolution(16)
	arrow.SetTipLength(0.3)
	arrow.SetTipRadius(0.1)

	glyph = vtk.vtkGlyph3D()
	glyph.SetSourceConnection(arrow.GetOutputPort())
	glyph.SetInputConnection(maskPts.GetOutputPort())
	glyph.SetVectorModeToUseNormal()
	glyph.SetScaleFactor(1)
	glyph.SetColorModeToColorByVector()
	glyph.SetScaleModeToScaleByVector()
	glyph.OrientOn()
	glyph.Update()
	return glyph

'''
	Make and display the surface.
	:param: st - the surface to display.
	:return The vtkRenderWindowInteractor.
'''
def DisplaySurface():

	# ------------------------------------------------------------
	# Create the surface, lookup tables, contour filter etc.
	# ------------------------------------------------------------
#	src = vtk.vtkPolyData()
	src = MakeParametricSource()
	# The scalars are named "Scalars"by default
	# in the parametric surfaces, so change the name.
	src.GetPointData().GetScalars().SetName("Elevation");
	scalarRange = src.GetScalarRange()

	lut = MakeLUT(total_bands)
	lut.SetTableRange(scalarRange)
	numberOfBands = lut.GetNumberOfTableValues()
	# bands = MakeIntegralBands(scalarRange)
	bands = MakeBands(scalarRange, numberOfBands, False)

	# Let's do a frequency table.
	# The number of scalars in each band.
	#print Frequencies(bands, src)

	# We will use the midpoint of the band as the label.
	labels = []
	for i in range(len(bands)):
		labels.append('{:4.2f}'.format(bands[i][1]))

	# Annotate
	values = vtk.vtkVariantArray()
	for i in range(len(labels)):
		values.InsertNextValue(vtk.vtkVariant(labels[i]))
	for i in range(values.GetNumberOfTuples()):
		lut.SetAnnotation(i, values.GetValue(i).ToString());

	# Create a lookup table with the colors reversed.
	# lutr = ReverseLUT(lut)

	# Create the contour bands.
	bcf = vtk.vtkBandedPolyDataContourFilter()
	bcf.SetInputData(src)
	# Use either the minimum or maximum value for each band.
	for i in range(0, numberOfBands):
		bcf.SetValue(i, bands[i][2])
	# We will use an indexed lookup table.
	bcf.SetScalarModeToIndex()
	bcf.GenerateContourEdgesOn()

	# Generate the glyphs on the original surface.
	#glyph = MakeGlyphs(src,False)

	# ------------------------------------------------------------
	# Create the mappers and actors
	# ------------------------------------------------------------
	srcMapper = vtk.vtkPolyDataMapper()
	srcMapper.SetInputConnection(bcf.GetOutputPort())
	srcMapper.SetScalarRange(scalarRange)
	srcMapper.SetLookupTable(lut)
	srcMapper.SetScalarModeToUseCellData()
 
	srcActor = vtk.vtkActor()
	srcActor.SetMapper(srcMapper)
	srcActor.RotateX(-45)
	srcActor.RotateZ(45)

	# Create contour edges
	edgeMapper = vtk.vtkPolyDataMapper()
	edgeMapper.SetInputData(bcf.GetContourEdgesOutput())
	edgeMapper.SetResolveCoincidentTopologyToPolygonOffset()

	edgeActor = vtk.vtkActor()
	edgeActor.SetMapper(edgeMapper)
	edgeActor.GetProperty().SetColor(0, 0, 0)
	edgeActor.RotateX(-45)
	edgeActor.RotateZ(45)
 
	glyphMapper = vtk.vtkPolyDataMapper()
	#glyphMapper.SetInputConnection(glyph.GetOutputPort())
	glyphMapper.SetScalarModeToUsePointFieldData()
	glyphMapper.SetColorModeToMapScalars()
	glyphMapper.ScalarVisibilityOn()
	glyphMapper.SelectColorArray('Elevation')
	# Colour by scalars.
	glyphMapper.SetScalarRange(scalarRange)
 
	glyphActor = vtk.vtkActor()
	glyphActor.SetMapper(glyphMapper)
	glyphActor.RotateX(-45)
	glyphActor.RotateZ(45)
 
	# Add a scalar bar.
	scalarBar = vtk.vtkScalarBarActor()
	# scalarBar.SetLookupTable(lut)
	# Use this LUT if you want the highest value at the top.
	scalarBar.SetLookupTable(lut) #lutr not lut to reverse
	scalarBar.SetTitle('Elevation (m)')
 
	# create a text actor
	txt = vtk.vtkTextActor()
	txt.SetInput(PrintColorBlindType())
	txtprop=txt.GetTextProperty()
	txtprop.SetFontFamilyToArial()
	txtprop.SetFontSize(18)
	txtprop.SetColor(0,0,0)
	txt.SetDisplayPosition(20,30)
 
	# ------------------------------------------------------------
	# Create the RenderWindow, Renderer and Interactor
	# ------------------------------------------------------------
	ren = vtk.vtkRenderer()
	renWin = vtk.vtkRenderWindow()
	renWin.SetWindowName("Test")
 	renWin.AddRenderer(ren)
	iren = vtk.vtkRenderWindowInteractor()
	iren.SetInteractorStyle(None)
	iren.SetRenderWindow(renWin)
 
	# add actors
	ren.AddViewProp(txt)
	ren.AddViewProp(srcActor)
	ren.AddViewProp(edgeActor)
	ren.AddViewProp(glyphActor)
	ren.AddActor2D(scalarBar)
 
	ren.SetBackground(0.7, 0.8, 1.0)
	renWin.SetSize(800, 800)
	renWin.Render()
 
	ren.GetActiveCamera().Zoom(1.5)
 
	return (ren, renWin, iren)


def PrintColorBlindType():
    if isColorblind == 0:
        return "NONE"
    elif isColorblind == 1:
        return "ALL"
    elif isColorblind == 2:
        return "DEUTERANOPIA"
    elif isColorblind == 3:
        return "PROTONOPIA"
    elif isColorblind == 4:
        return "TRITANOPIA"

#    return {
#        "NONE": isColorblind == 0,
#        "ALL": isColorblind == 1,
#        "DEUTERANOPIA": isColorblind == 2,
#        "PROTONOPIA": isColorblind == 3,
#        "TRITANOPIA": isColorblind == 4
#    }


''' 

	Custom interactions 

'''


'''
	Handle the mouse button events.
'''
def ButtonEvent(obj, event):
    global Rotating, Panning, Zooming
    if event == "LeftButtonPressEvent":
        Rotating = 1
    elif event == "LeftButtonReleaseEvent":
        Rotating = 0
    elif event == "MiddleButtonPressEvent":
        Panning = 1
    elif event == "MiddleButtonReleaseEvent":
        Panning = 0
    elif event == "RightButtonPressEvent":
        Zooming = 1
    elif event == "RightButtonReleaseEvent":
        Zooming = 0

'''
	General high-level logic
'''
def MouseMove(obj, event):
    global Rotating, Panning, Zooming
    global iren, renWin, ren
    lastXYpos = iren.GetLastEventPosition()
    lastX = lastXYpos[0]
    lastY = lastXYpos[1]

    xypos = iren.GetEventPosition()
    x = xypos[0]
    y = xypos[1]

    center = renWin.GetSize()
    centerX = center[0]/2.0
    centerY = center[1]/2.0

    if Rotating:
        Rotate(ren, ren.GetActiveCamera(), x, y, lastX, lastY,
               centerX, centerY)
    elif Panning:
        Pan(ren, ren.GetActiveCamera(), x, y, lastX, lastY, centerX,
            centerY)
    elif Zooming:
        Dolly(ren, ren.GetActiveCamera(), x, y, lastX, lastY,
              centerX, centerY)

'''
	The user pressed a key, so we need to handle
	the event appropriately.
'''
def Keypress(obj, event):
	global isColorblind, DEUTERANOPIA, PROTONOPIA, TRITANOPIA

	key = obj.GetKeySym()
	if key == "e":
		obj.InvokeEvent("DeleteAllObjects")
		sys.exit()
	elif key == "w":
		Wireframe()
	elif key =="s":
		Surface()
	elif key =="a":
		isColorblind = ALL
		UpdateColor()
	elif key =="d":
		isColorblind = DEUTERANOPIA
		UpdateColor()
	elif key =="p":
		isColorblind = PROTONOPIA
		UpdateColor()
	elif key =="t":
		isColorblind = TRITANOPIA
		UpdateColor()
	elif key =="space":
		isColorblind = 0
		UpdateColor()

'''
	If the user pressed a key corresponding
	to a change in the colorblind options
	we need to redraw the scene.
	
	NOTE: acts like glutredisplay();
'''
def UpdateColor():
    pick_colors()
    lut = MakeLUT(total_bands)
    renWin = iren.GetRenderWindow()
    rens = renWin.GetRenderers()
    rens.InitTraversal()
    ren_temp = rens.GetNextItem()
    while ren_temp:
        actors = ren_temp.GetActors()
        actors.InitTraversal()
        actor = actors.GetNextItem()
        while actor:                
            actor.GetMapper().SetLookupTable(lut)
            actor = actors.GetNextItem()
        actors2D = ren_temp.GetActors2D()
        actors2D.InitTraversal()
        actor2D = actors2D.GetNextItem()
        while actor2D:
            if actor2D.IsA("vtkTexturedActor2D") == True:
                actor2D.SetInput(PrintColorBlindType())
            else:    
                actor2D.SetLookupTable(lut)
            actor2D = actors2D.GetNextItem()
        ren_temp = rens.GetNextItem()
    renWin.Render()
    
'''
	Routines that translate the events into camera motions.

	This one is associated with the left mouse button. It translates x
	and y relative motions into camera azimuth and elevation commands.
'''
def Rotate(renderer, camera, x, y, lastX, lastY, centerX, centerY):
    camera.Azimuth(lastX-x)
    camera.Elevation(lastY-y)
    camera.OrthogonalizeViewUp()
    renWin.Render()


'''
	Pan translates x-y motion into translation 
	of the focal point and position.
'''
def Pan(renderer, camera, x, y, lastX, lastY, centerX, centerY):
    FPoint = camera.GetFocalPoint()
    FPoint0 = FPoint[0]
    FPoint1 = FPoint[1]
    FPoint2 = FPoint[2]

    PPoint = camera.GetPosition()
    PPoint0 = PPoint[0]
    PPoint1 = PPoint[1]
    PPoint2 = PPoint[2]

    renderer.SetWorldPoint(FPoint0, FPoint1, FPoint2, 1.0)
    renderer.WorldToDisplay()
    DPoint = renderer.GetDisplayPoint()
    focalDepth = DPoint[2]

    APoint0 = centerX+(x-lastX)
    APoint1 = centerY+(y-lastY)

    renderer.SetDisplayPoint(APoint0, APoint1, focalDepth)
    renderer.DisplayToWorld()
    RPoint = renderer.GetWorldPoint()
    RPoint0 = RPoint[0]
    RPoint1 = RPoint[1]
    RPoint2 = RPoint[2]
    RPoint3 = RPoint[3]

    if RPoint3 != 0.0:
        RPoint0 = RPoint0/RPoint3
        RPoint1 = RPoint1/RPoint3
        RPoint2 = RPoint2/RPoint3

    camera.SetFocalPoint( (FPoint0-RPoint0)/2.0 + FPoint0,
                          (FPoint1-RPoint1)/2.0 + FPoint1,
                          (FPoint2-RPoint2)/2.0 + FPoint2)
    camera.SetPosition( (FPoint0-RPoint0)/2.0 + PPoint0,
                        (FPoint1-RPoint1)/2.0 + PPoint1,
                        (FPoint2-RPoint2)/2.0 + PPoint2)
    renWin.Render()


# Dolly converts y-motion into a camera dolly commands.
def Dolly(renderer, camera, x, y, lastX, lastY, centerX, centerY):
    dollyFactor = pow(1.02,(0.5*(y-lastY)))
    if camera.GetParallelProjection():
        parallelScale = camera.GetParallelScale()*dollyFactor
        camera.SetParallelScale(parallelScale)
    else:
        camera.Dolly(dollyFactor)
        renderer.ResetCameraClippingRange()

    renWin.Render()

'''
	Wireframe sets the representation of all actors to wireframe.
'''
def Wireframe():
    actors = ren.GetActors()
    actors.InitTraversal()
    actor = actors.GetNextItem()
    while actor:
        actor.GetProperty().SetRepresentationToWireframe()
        actor = actors.GetNextItem()

    renWin.Render()

'''
	Surface sets the representation of all actors to surface.
'''
def Surface():
    actors = ren.GetActors()
    actors.InitTraversal()
    actor = actors.GetNextItem()
    while actor:
        actor.GetProperty().SetRepresentationToSurface()
        actor = actors.GetNextItem()
    renWin.Render()	

'''
	Randomly choose the rgb values.
'''
def choose_rgb(val):

	global dispColor1, dispColor2, isColorblind	
	temp = Color(0.0,0.0,0.0) # temporary color
	
	notsat = 0
	if isColorblind == DEUTERANOPIA:
		notsat = 1 # green
	elif isColorblind == PROTONOPIA:
		notsat = 0 # red
	elif isColorblind == TRITANOPIA:
		notsat = 2 # blue
	else:
		notsat = randint(0,2) # randomly pick the value of 0 for rgb
	
	sat = randint(0,2) # randomly pick the value of 1 for rgb
	other = randint(0,2)
	while notsat == sat or notsat == other or other == sat:
		if notsat == sat:
			sat = randint(0,2)
		if notsat == other:
			other = randint(0,2)
		if other == sat:
			other = randint(0,2)
	
	# set rgb of temp Color
	temp.set_value(sat, 1.0)
	temp.set_value(notsat, 0.0)
	fval = random.uniform(0.0, 1.0)
	
	# try and get rid of the yellow color
	if isColorblind == TRITANOPIA and other == 1:
		while fval <= 0.7:
			fval = random.uniform(0.0, 1.0)
	
	temp.set_value(other, fval)
	
	# set the first color
	if val == 1:
		dispColor1.set_red(temp.r)
		dispColor1.set_green(temp.g)
		dispColor1.set_blue(temp.b)
	# set the second color
	elif val == 2:
		dispColor2.set_red(temp.r)
		dispColor2.set_green(temp.g)
		dispColor2.set_blue(temp.b)

	
'''
	Randomly pick two colors. Ensure that the colors
	are both highly saturated and that they are a certain
	distance apart.
''' 
def pick_colors():

	global dispColor1, dispColor2, isColorblind
	
	if isColorblind == ALL:
		# black color
		dispColor1.set_red(0.0)
		dispColor1.set_green(0.0)
		dispColor1.set_blue(0.0)			
		# white color
		dispColor2.set_red(1.0)
		dispColor2.set_green(1.0)
		dispColor2.set_blue(1.0)
	else:
		choose_rgb(1)
		choose_rgb(2)
	
		times_compared = 0
		dist = dispColor1.get_distance(dispColor2)
		while dist <= 1.0:
			#update the number of times we have compared colors
			times_compared += 1
		
			if times_compared >= 10:
				# pick a new first color
				choose_rgb(1)
				
				# reset the number of times compared
				times_compared = 0
				
			# pick a new second color
			choose_rgb(2)
			
			# recalculate the distance
			dist = dispColor1.get_distance(dispColor2)
	
'''
	Create a list of all of the colors that 
	will be displayed on the range. This is 
	one of two algorithms that will be used to 
	display the colors.
'''
def select_colors(num_distinct):
	the_colors = []
	the_colors.append([0.0, dispColor1.r, dispColor1.g, dispColor1.b])
	
	r_dist = abs(dispColor1.r - dispColor2.r) / float(num_distinct - 2)
	g_dist = abs(dispColor1.g - dispColor2.g) / float(num_distinct - 2)
	b_dist = abs(dispColor1.b - dispColor2.b) / float(num_distinct - 2)
 
	i = 0
	while i < (num_distinct-2):
		red = 0.0
		green = 0.0
		blue = 0.0
		
		val = float(i + 1) / float(num_distinct)
		
		if dispColor1.r > dispColor2.r:
			red = dispColor1.r - r_dist
		else:
			red = dispColor1.r + r_dist
			
		if dispColor1.g > dispColor2.g:
			green = dispColor1.g - g_dist
		else:
			green = dispColor1.g + g_dist
		
		if dispColor1.b > dispColor2.b:
			blue = dispColor1.b - b_dist
		else:
			blue = dispColor1.b + b_dist
		
		the_colors.append([val, red, green, blue])
		i += 1
	
	
	the_colors.append([1.0, dispColor2.r, dispColor2.g, dispColor2.b])
	
	return the_colors
 
#if __name__ == '__main__':
	'''
	Determine what arguments were provided by the user.
	Also provide the user with some information about the 
	program using the -help option.
	'''
if len(sys.argv) > 1:
	if sys.argv[1] == "-h" or\
		sys.argv[1] == "--h" or\
		sys.argv[1] == "--help" or \
		sys.argv[1] == "-help":
		print("\n-cb 	\tColorblind option.")
		print("\tall : \tCompletely Colorblind")
		print("\tnone : \tNot Colorblind")
		print("\td : \tDeuteranopia (Colorblind in the Medium Wavelength)")
		print("\tp : \tProtonopia   (Colorblind in the Long   Wavelength)")
		print("\tt : \tTritanopia   (Colorblind in the Short  Wavelength)")
		print("\n-b 	\tNumber of bands option. Controls the number of distinct bands will be displayed.")
		sys.exit(1)
	else:
		if len(sys.argv) != 1 and len(sys.argv) != 3 and len(sys.argv) != 5:
			print("Please provide the correct arguments and options.")
			print("Use the help option to determine appropriate options.")
			sys.exit(1)
	
		i = 1 # start at the first optional arguement
		while i < len(sys.argv):
			if sys.argv[i] == "-cb" and sys.argv[i+1] == "all":
				# completely colorblind
				isColorblind = ALL
			elif sys.argv[i] == "-cb" and sys.argv[i+1] == "d":
				# color blind to the medium wavelength
				isColorblind = DEUTERANOPIA
			elif sys.argv[i] == "-cb" and sys.argv[i+1] == "p":
				# color blind to the long wavelength
				isColorblind = PROTONOPIA
			elif sys.argv[i] == "-cb" and sys.argv[i+1] == "t":
				# color blind to the short wavelength
				isColorblind = TRITANOPIA
			elif sys.argv[i] == "-cb" and sys.argv[i+1] == "none":
				# color blind to the short wavelength
				isColorblind = 0
				pick_colors()
			elif sys.argv[i] == "-cb":
				print("Invalid colorblind option: Using NONE.")
			elif sys.argv[i] == "-b":
				if int(sys.argv[i+1]) >= 3:
					total_bands = int(sys.argv[i+1])
				else:
					print("Invalid number of bands: Using 8.")
			i = i + 2
pick_colors()

 # Add the observers to watch for particular events. These invoke
# Python functions.
Rotating = 0
Panning = 0
Zooming = 0 
 
ren, renWin, iren = DisplaySurface()
 
iren.AddObserver("LeftButtonPressEvent", ButtonEvent)
iren.AddObserver("LeftButtonReleaseEvent", ButtonEvent)
iren.AddObserver("MiddleButtonPressEvent", ButtonEvent)
iren.AddObserver("MiddleButtonReleaseEvent", ButtonEvent)
iren.AddObserver("RightButtonPressEvent", ButtonEvent)
iren.AddObserver("RightButtonReleaseEvent", ButtonEvent)
iren.AddObserver("MouseMoveEvent", MouseMove)
iren.AddObserver("KeyPressEvent", Keypress)
	
	
iren.Initialize()
renWin.Render()
iren.Start()
 
#	my_iren.Render()
#	
#	my_iren.AddObserver("KeyPressEvent", Keypress)
#	
#	my_iren.Start()
#	my_iren.Initialize()
