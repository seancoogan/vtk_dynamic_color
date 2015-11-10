import math
import sys
import os

''' 
The color class will contain the red, green, and blue
values that are associated with the color. These rgb
values must be between 0 and 1.
 '''
class Color:
	r = 0.0 # red portion
	g = 0.0 # green portion
	b = 0.0 # blue portion
	
	brightness = 0.0 # brightness dependent on all rgb values
	
	# create a constructor with some default values
	# in case the user never passes the values
	def __init__(self, r=0.0, g=0.0, b=0.0):
		# make sure that the red value is a float
		if isinstance(r, float) and r <= 1.0 and r >= 0.0:
			self.r = r # set the red portion to r
		else:
			print("Bad RED Value: Red = 0.0")
			self.r = 0.0
			
		# make sure that the green value is a float
		if isinstance(g, float) and g <= 1.0 and g >= 0.0:
			self.g = g # set the green portion to g
		else:
			print("Bad GREEN Value: Green = 0.0")
			self.g = 0.0
			
		# make sure that the blue value is a float
		if isinstance(b, float) and b <= 1.0 and b >= 0.0:
			self.b = b # set the blue portion to b
		else:
			print("Bad BLUE Value: Blue = 0.0")
			self.b = 0.0
	
	# set the red parameter of the color
	def set_red(self, r):
		# make sure that the red value is a float
		if isinstance(r, float) and r <= 1.0 and r >= 0.0:
			self.r = r # set the red portion to r
		else:
			print("Bad RED Value: Red = 0.0")
			self.r = 0.0

	# set the green parameter of the color
	def set_green(self, g):
		# make sure that the green value is a float
		if isinstance(g, float) and g <= 1.0 and g >= 0.0:
			self.g = g # set the green portion to g
		else:
			print("Bad GREEN Value: Green = 0.0")
			self.g = 0.0
	
	# set the blue parameter of the color
	def set_blue(self, b):
		# make sure that the blue value is a float
		if isinstance(b, float) and b <= 1.0 and b >= 0.0:
			self.b = b # set the blue portion to b
		else:
			print("Bad BLUE Value: Blue = 0.0")
			self.b = 0.0
	
	# set the arbitrary value in the Color
	def set_value(self, index, value):
		if not isinstance(value, float) or \
			value > 1.0 or \
			value < 0.0:
			print("Check the index or value!")
			return
	
		if index == 0:
			self.r = value
		elif index == 1:
			self.g = value
		elif index == 2:
			self.b = value
		else:
			print("Bad index!")
			return
	
	# calculate the distance between this color
	# and another color
	def get_distance(self, other):
		# initial error checking
		if not isinstance(other, Color):
			# we should return -- Error
			print("Please only compare a color.")
			return
		
		# calculate the differences between the colors
		difr = self.r - other.r
		difg = self.g - other.g
		difb = self.b - other.b
		# return the distance formula
		return math.sqrt(math.pow(difr,2) + math.pow(difg,2) + math.pow(difb,2))
		
	# calculate the brightness of this color
	# based on all RGB values
	#
	# Based on a weighted sum
	def calculate_bightness(self):
		self.brightness = math.sqrt( \
					(.299 * math.pow(self.r,2)) + \
					(.587 * math.pow(self.g,2)) + \
					(.114 * math.pow(self.b,2)))
	
		return self.brightness
	
	# calculate the midpoint color coordinates between 
	# this color and the color that is passed in
	def calculate_midpoint(self, other):#, num_distinct):
		# initial error checking
		if not isinstance(other, Color):
			# we should return -- Error
			print("Please only compare a color.")
			return
			
		# calculate that mid point
		mid_r = ((self.r + other.r) / 2) # + 1/num_distinct
		mid_g = ((self.g + other.g) / 2) #+ 1/num_distinct
		mid_b = ((self.b + other.b) / 2) #+ 1/num_distinct
		
		return (mid_r, mid_g, mid_b)
		
