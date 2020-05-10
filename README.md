# Simple Cosmic String Network Simulation
Cosmic strings are topological defects that likely form during phase transitions in the early universe. In flat space-times, these strings obey simple dynamics--little more than the ordinary 1D wave-equation. This program considers some number of strings generated inside a cubic (spacial) volume, and animates their motion over (conformal) time. The algorthim is adapted from from Smith & Vilenkin (Physical Review D, Vol 36, Number 4), albeit with more simplistic (and unrealistic) boundary conditions and a non-physical string generation walk.

## Features
Allows for one to generate arbitrary number of strings in a fixed box, and animate their motion over time. MP4 files can be saved. Graphical/presentational tweaks can be made to suit one's needs.

## Dependencies
The script was written using Python 3.7, and imports from Numpy and Matplotlib libraries. 


## Parameters
User-defined variables include:
- Number length of box `N` (int)
- Space-time step size `delta`(float)
- How many times to repeat a step in the generating walk `rep` (int)
- Whether or not to run self-intersection routine `self_intersect` (bool)
- Whether or not to save an .mp4 video `save_anim` (bool)
- Number of strings in the network `number_strings` (int)
- Number of frames to simulate (number of time-steps) `number_frames` (int)
- Framerate `framerate`(float) of output .mp4

## Classes
There are two central objects in the program, belonging to either the `CosmicString` class or `CosmicStringLoop` class. `CosmicString` initializes by taking a random walk to generate the spacial coordinates of a string segment. It then calls `dynamics()` to calculate the motion of each coordinate. If `self_intersect = True`, then `dynamics()` in `CosmicString` will enter a routine to check for overlapping coordinates, and initialize objects in the `CosmicStringLoop` class. Everything ultimately follows the same dynamics. 

*Warning: Enabling self-intersection in its current state is expensive.*

## Liscense
This project is filed under the MIT liscense.
