NPVol
=====
This program is designed to find the volume of a nanoparticle using Monte Carlo
integration. The volume of your nanoparticle is a nice thing to know in certain
cases, because you can use it to

 * apply a size-correction factor to your dielectric constant to account for
   quantum-size effects
 * determine the absorption/scattering/extinction cross-section of your nanoparticle
 * determine the density of your nanoparticle
 
I'm sure there are many other wonderful uses, but these are why I am interested
in the volume. Please see the **Theory** section at the bottom for an
explaination of how this code works.

Compiling NPVol
---------------

Compiling `NPVol` should be as easy as running `Make` in the `npvol` directory.
If you wish to install the binary to your PATH, then you can also run `make
install`.  There are a few customizations that I will go over for you.

I have deliberately avoided a configure script or any other heavy duty build
system in favor of a simple Makefile for this project.  The reason is that
there are no external dependencies (except for MPI if you choose to use it) so
there is no need to use excessive force.  As such, I have added in some
customizations to the `Makefile` itself to simulate on a small scale how these
configure systems work.  Feel free to change what needs changing, or better yet
improve and generalize what is there.

### Known Compilers

The `Makefile` for `NPVol` is aware of the compilers `gfortran`, `ifort`,
`g95`, and `pgf90`, and thus has compilation flags for these compilers
pre-programmed.  The default compiler is `gfortran` (since this is what most
people have), but I have found the `ifort` gives far superior performance, so
use that if you have access (`g95` is particulary bad, so I don't recommed).
To choose the compiler you want, use

	make FC=ifort

where in this case it is assumed you wanted `ifort`.

### Parallelization Methods

By default a serial version of `NPVol` is compiled.  This is fine if your
nanoparticles are less than about 2000 atoms, but for particles larger than
this you will want to use parallelized code.  `NPVol` has been written to
support either **OpenMP** or **MPI**, but not both.  

##### Which method should I choose?

In all cases, I recommend **OpenMP**.  It is built into the compiler, so you
don't need to install special libraries to use it.  Also, you execute the
program the same way as the serial version.  

The only time I recommend **MPI** is if you are on a distributed memory system
such as a supercomputer or using `g95`, since neither will support **OpenMP**
(**OpenMP** is shared memory only, and `g95` hasn't implemented it yet).

##### Compiling a parallel method

To compile with parallelization support, you must set the `PARALLEL` variable.
For **OpenMP**:

	make PARALLEL=openmp

For **MPI**:

	make PARALLEL=mpi

**MPI** comes in many *flavors* (implementations), such as **OpenMPI** (not to
be confused with **OpenMP**), **MPICH**, etc.  `NPVol`'s `Makefile` is aware of
**OpenMPI** and **MPICH**, as these are probably the two most common free
**MPI** flavors.  You should determine which flavor you have installed on your
machine before you compile `NPVol`.  By default, the `Makefile` assumes you
have **OpenMPI** installed.  To change this, use the `FLAVOR` variable

	make PARALLEL=mpi FLAVOR=mpich

If you have a different flavor, please add it to the `Makefile`!

Installing NPVol
----------------

If you don't want your binary to be in the extracted folder, you may `install`
the binary to your hard drive:

	make install
	
The default `PREFIX` value is `/usr/local/`, so by default this will install
the `NPVol` binary to `/usr/local/bin`.  If you don't have write permissions to
this directory, you may have to use `sudo`:

	sudo make install
	
If you do not have `sudo` access, the I recommend installing to the `bin`
directory in your `$HOME` folder:

	make install PREFIX=$HOME

If `$HOME/bin` (`~/bin`) is not in your `$PATH`, then you can add it by adding
this line to your `.bashrc` (or `.profile` on Mac OS X):

	export PATH=$PATH:$HOME/bin

Using NPVol
-----------

To execute `NPVol`, you simply must give it an .xyz file as input

	NPVol nanoparticle.xyz
	
where it is assumed that you have `install`ed `NPVol` to your `$PATH`.  

### What is the format of an XYZ file?

An XYZ file simply gives the atomic name and coordinates of each atom in your
system, with the number of atoms total on the first line and a comment on the
second:

	3
	This is a comment
	Au  0.0 0.0 0.0
	Ag  3.0 0.0 0.0
	Ag -3.0 0.0 0.0

Note that the coordinates are Cartesian coordinates and by default are assumed
to be in Angstroms, although this may be overridden.  For `NPVol` an
additional restriction is added in that the comment line must give the radii of
each atom type in the system:

	3
	Ag=1.4445 Au=1.4445
	Au  0.0 0.0 0.0
	Ag  3.0 0.0 0.0
	Ag -3.0 0.0 0.0
	
Note that you must have an equals between the atom name and the radius and
there can be no spaces.

### Options

`NPVol` has a few options

 * `-h` Displays a help message and quits
 * `-?` Displays a help message and quits
 * `-b` Reads in coordinates and radii in Bohr instead of Angstroms
 * `-s NSAMPLES` Override the default number of sample points to a chosen
    value.  The default number of samples is given by
     **# Samples = MAX(# Atoms * 50, 10^6)**.
 * `-o` Optimize the radius multiplier.  This is explained below.
 * `-t THRESHOLD` Threshold for optimization.  This is explained below.

##### Radius multiplier optimization

Because the number of sample points must increase as the system size grows, the
timing of `NPVol` scales poorly with system size (in fact is N^2). This is fine
for a stand-alone code such as `NPVol`, since the idea is that you may use
`NPVol` to calculate the volume of your nanoparticle and then use that value as
input in other codes that are aimed at calculating other properties but need
the volume.  

For most systems this poor scaling is barely noticable, since the overhead to
`NPVol` is quite low.  The problem is that when the number of atoms becomes
very large then `NPVol` takes a very long time.  For example, for 1070000
atoms, it takes `NPVol` 8 hours to caluculate the volume using the default
settings and running on 8 processors; that would be 64 hours if run in serial!
You can always use more processors ( `NPVol` is *embarassingly parallel*,
meaning that it scales perfectly with # processors) but not everybody has
access to a large number of processors.

To circumvent this, some codes that need the volume give the option of simply
summing up the volumes of the atoms assuming they are spheres with the given
radii.  This method is **_much_** faster than Monte Carlo, but it will
underestimate the volume since it does not take into account the space between
atoms (which is accounted for in `NPVol`).  To compensate, the radii may be
multiplied by some factor.  Using the `-o` option will cause `NPVol` to
calculate the optimal multiplier that you may give so that any nanoparticle's
volume with similar shape, structure and composition can be calculated quickly.
This optimization is done by actually calculating the volume using the
multiplier and then determining the error between the two methods.  The default
threshold for the error is 1%, but this can be edited using the `-t` option.

### Running in parallel

##### OpenMP

As mentioned previously, if you compiled with **OpenMP** then there is no need
to execute this code any differently than you would if it were serial:

	NPVol nanoparticle.xyz

However, if you wish to use less processors than the maximum on your machine,
you may set the `OMP_NUM_THREADS` environment variable.  Let's say you only
wanted to run on 2 processors:

	OMP_NUM_THREADS=2 NPVol nanoparticle.xyz
	
##### MPI

Running with **MPI** requires that you use **MPI**'s executable to set up the
parallel environment.  Let's say you wanted to run on 4 processors:

	mpirun -n4 NPVol nanoparticle.xyz
	
You **must** give the number of processors.  In somecases, `mpirun` is not
available, but `mpiexec` is an equivalent routine:

	mpiexec -n4 NPVol nanoparticle.xyz

Theory
------

Imagine that you want to know the value of pi, and all you have to determine
this is a circular dartboard inside a square frame and an infinite number of
darts.  It turns out that this is all we need to determine pi (honest).  

We know that the area of the square frame that the dartboard is in is 

  A^square = (2r)^2 = 4r^2

where r is half the side length of the frame (you'll see why it was defined
that way in a moment). We also know that the area of the circular dartboard is 

  A^circle = &pi;r^2

Now, if there was some way to find the ratio (R) of the areas, then we could
find pi, like this:

  R = A^circle / A^square = ( &pi;r^2 ) / ( 4r^2 ) = &pi; / 4  
  &pi; = 4R

This means that &pi; is just 4 times the ratio of the two areas!

…But wait a minute.  What does that have to do with my dartboard?  How can I
find the area of the circle without pi?  This is where the Monte Carlo
technique comes in. In Monte Carlo, we use random numbers to compare something
known (or easily determined) to something unknown (or less easily determined).  

Remember that our circular dartboard is enclosed by a square frame, with the
side of the frame being exactly equal to the diameter of the dartboard.  If you
are a terrible dart player, and every dart you throw has a random chance of
hitting *anywhere* in the frame, be it on the dartboard or not, then the ratio
of darts that are in the frame to darts that are only within the dartboard is
directly proportional to the ratio of the areas, and therefore is equal to the
R defined above!.  

I'll let that soak in a bit. Try drawing a circle enclosed by a square on a
piece of paper, then putting random dots *in the square, as if the circle were
not there*.  This is exactly what I am talking about.  The more points you put
(darts you throw), the closer your ratio is to exactly the ratio of the areas.
If you were to use an infinite number of points/darts, then you would get the
exact answer.

So, how does this relate to nanoparticles?  In this case, the dartboard is the
nanoparticle, and the frame is a box that completely surrounds the
nanoparticle. The volume of the box is easy to calculate, but obviously the
volume of an arbitrarily shaped nanoparticle is not.  All we know about our
nanoparticle is the atomic positions in space and the radii of each of these
atoms.  If we take random sample points inside the box, and take the ratio of
sample points that were in the box to those only in the nanoparticle, then we
can find the volume of the nanoparticle:

  R = points^np / points^box  
  V^np = R * V^box
	
Thus, we can easily find the volume of our arbitrarily shaped nanoparticle!
