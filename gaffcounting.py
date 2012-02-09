import argparse
import sys
from statemodel import *
from statesolver import StateSolver


'''The first section just parses command-line arguments.'''

parser = argparse.ArgumentParser(description='Finds Gaffnian states for given n, and N. Counts number of states and compares to the counting forumla.')

parser.add_argument('n', type=int, help='Specify domain wall number.')
parser.add_argument('N', type=int, help='Specify particle number.')
parser.add_argument('-s', action='store_true', default=False, \
					help='Sphere option, default False. '+\
					'If this flag is set, will find states on the sphere and not on torus.')
parser.add_argument('-v', action='store_true', default=False, \
					help='Verbose, default False. If set, all patterns will be printed.')
parser.add_argument('-z', action='store_true', default=False, \
					help='Debugging (extra verbose), default False. '+\
					'If set, statements will be printed at all steps of state-checking recursion.')


args = parser.parse_args(sys.argv[1:])

n = args.n
N = args.N
sphere_flag = args.s
verbose_flag = args.v
debug_flag = args.z

if n%2 != N%2:
	sys.exit('Please choose n and N either both even or both odd. Run with -h for help.')


'''Creates a Gaffnian model with the given n,N on the proper geometry.
If you want to write your own model, put its name here.'''
if sphere_flag:
	state_model = GaffnianSphere(n,N)
else:
	state_model = GaffnianTorus(n,N)



'''Creates the StateSolver, which will find all the states and store them.
This should be general, as long as your StateModel has the right information.'''
solver = StateSolver( state_model, debug_flag )
solved_states = solver.state_list

if verbose_flag:
	print 'States on '+state_model.geometry
	print solver.print_states(solved_states)

print 'Number of states on {}: {}'.format(state_model.geometry,len(solved_states))


'''Finds the output of the counting formula (supplied in the StateModel).'''
if verbose_flag: print 'L = {}, Fmax = {}, Fmin = {}'.format(state_model.L, state_model.Fmax, state_model.Fmin)

predicted_number_of_states = 0
for F in range(state_model.Fmax,state_model.Fmin -1,-2):
	if verbose_flag: 
		print 'F = {} states from {} formula:'.format(F, state_model.geometry)
		print '{}*{} = {}*{} = {}'.format( \
			state_model.pos_degen_string(F),state_model.internal_degen_string(F),\
			state_model.pos_degen(F),state_model.internal_degen(F),state_model.count(F))
	predicted_number_of_states += state_model.count(F)

print 'Number of states from {} formula: {}'.format(state_model.geometry, predicted_number_of_states)
