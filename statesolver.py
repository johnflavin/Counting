class StateSolver:
	def __init__(self, model, debug=False):
		self.L = model.L						# number of orbitals
		self.k = model.cell_size				# at most max_in_k particles in k sites
		self.N = model.N						# number of particles
		self.max_in_k = model.max_contiguous	# at most max_in_k particles in k sites
		self.max = model.max_occupation			# at most max_occupation particles in any one site
		self.geometry = model.geometry			# string holding the geometry (e.g. 'torus', 'sphere')
		self.debug = debug						# prints messages throughout recursion if set
		
		self.state_list = self.find_states( State(self.L,self.k,self.geometry) )

	def find_states(self, state):
		if self.debug: print 'find_states called with {} focused on index {}'.format(state.string(),state.focus)
		if state.running_sum==self.N:
			if self.debug: print 'state is good. passing it up.'
			return [state]
		if state.focus==self.L:
			if self.debug: print 'reached the end and state is still not good'
			return []

		children = []
		for val in range(self.max+1):
			if self.debug: print 'considering {} in position {} in {}...'.format(val,state.focus,state.string())
			prox = state.proximate_list()
			if self.debug: print 'proximity numbers around index {}: {}'.format(state.focus,prox)
			new_prox = [p+val for p in prox]
			if self.debug: print 'if we add value {} at index {}, prox. numbers become {}'.format(val,state.focus,new_prox)
			
			if max(new_prox) <= self.max_in_k:
				if self.debug: print 'ok. adding {} to {} at index {}'.format(val,state.string(),state.focus)
				newstate = state.clone()
				newstate.change_value(val)
				if self.debug: print 'changing focus to {}'.format(state.focus +1)
				newstate.change_focus(state.focus+1)

				if self.debug: print 'passing down '+newstate.string()
				output_from_below = self.find_states(newstate)
				if self.debug: 
					print 'level {} - received from below:'.format(state.focus)
					print self.print_states(output_from_below)

				if len(output_from_below)>0:
					children.extend(output_from_below)

				if self.debug: 
					print 'children now:'.format(state.focus)
					print self.print_states(children)
			else:
				if self.debug:
					print 'value {} at index {} would cause an overflow'.format(val,state.focus)
				break 

		if self.debug: 
			print 'end of level {}. passing up '.format(state.focus)
			print self.print_states(children)
		return children

	def print_states(self, input_array):
		screen_width = 100
		num_states = len(input_array)

		if num_states==0:
			return '[]'
		
		group_length = self.k
		number_of_groups = self.L/group_length
		number_left_over = self.L-number_of_groups*group_length
		state_length = 3+number_of_groups*(group_length+1)+number_left_over

		num_states_per_line = screen_width / state_length
		num_lines = num_states / num_states_per_line
		num_on_last_line = num_states - num_lines*num_states_per_line

		output = ''
		for line in range(num_lines):
			output+=' '.join( [state.string() for state in \
							input_array[num_states_per_line*line:num_states_per_line*(line+1)]] )
			if num_lines-line>1: output+='\n'
		if num_on_last_line > 0:
			if num_lines>0: output+='\n'
			output+=' '.join( [state.string() for state in input_array[-num_on_last_line:]] )

		return output
		

class State:
	def __init__(self, pattern_size, cell_size, geometry='torus'):
		self.L = pattern_size
		self.k = cell_size
		self.geometry = geometry				# holds a string, e.g. 'torus' or 'sphere'

		self.pattern = [0]*pattern_size			# length L array holding the occupation number pattern.
		
		'''Proximity Array
			A length L array holding the sum of the occupation numbers of the k orbitals 
			around the one indexed. For odd k it is symmetric about the indexed element, but for
			even k it defaults to holding one extra element to the left.
			Examples: If k = 2 and pattern = [0,1,2,3,4,5,6], then 
			proximity = [6,1,3,5,7,9,11] on the torus and [0,1,3,5,7,9,11] on the sphere.
			If k = 3 and pattern is as above, then
			proximity = [7,3,6,9,12,15,11] on the torus and [1,3,6,9,12,15,11] on the sphere.

			If you want to include another geometry, the only thing you need to change is 
			proximate_indices. 
		'''
		self.proximity = [0]*pattern_size

		self.focus = 0							# The index that is currently being worked on
		self.running_sum = 0					# Total number of occupied sites in the pattern

	def clone(self):
		newstate = State(self.L,self.k,self.geometry)
		newstate.pattern = list(self.pattern)
		newstate.proximity = list(self.proximity)
		newstate.focus = self.focus
		newstate.running_sum = self.running_sum

		return newstate

	def change_focus(self, new_focus):
		self.focus = new_focus

	def change_value(self, value):
		self.pattern[self.focus] = value
		self.running_sum += value

		for index in self.proximate_indices():
			self.proximity[index] += value

	def proximate_list(self):
		return [self.proximity[index] for index in self.proximate_indices()]

	def proximate_indices(self):
		if self.k % 2 == 1:
			left_limit = self.focus - (self.k-1)/2
			right_limit = self.focus + (self.k-1)/2
		else:
			left_limit = self.focus - self.k/2
			right_limit = self.focus + self.k/2 - 1

		if self.geometry=='torus': 
			if right_limit >= self.L:
				right_limit = right_limit - self.L
				return range(right_limit+1)+range(left_limit,self.L)
		elif self.geometry=='sphere':
			if left_limit < 0:
				left_limit = 0
			if right_limit >= self.L:
				right_limit = self.L - 1
		# if you include another geometry, write how indices must be handled at the edges
		# elif self.geometry == 'your geometry':
		#	blah blah blah

		return range(left_limit,right_limit + 1)

	def string(self):
		group_length = self.k
		number_of_groups = self.L/group_length
		number_left_over = self.L-number_of_groups*group_length
	
		grouped = []
		for g in range(number_of_groups):
			grouped.append(''.join([str(value) for value in self.pattern[g*group_length:(g+1)*group_length]]))
		if number_left_over > 0:
			grouped.append(''.join([str(value) for value in self.pattern[-number_left_over:]]))
		
		return '['+' '.join(grouped)+']'