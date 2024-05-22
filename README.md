# imfdFlap
**tl;dr:**  
g++ imfdFlap.cpp  
./a.out prepare  
./allrun.sh  

Build imfdFlap.cpp using your c++ compiler of choice (I used g++ 11.4.0). Run the output file with "prepare" and execute the allrun.sh script. This creates a couple of sub-directories prefaced with a . and one sub-directory *cases*, which contains all preCICE coupling directories and their respective OpenFoam sub-directories. Any changes made to .folders can be pushed to the case directories using the copyFiles.sh script.

**Copy of the default output message:**  
Self-installing OpenFOAM experiment: IMFD-Flap using preCICE, by Alrik Luda - rev. 0 - (2024-04-14)

**About:**  
This file generates scripts to install the latest versions of OpenFOAM, OpenFOAM-Adapter, preCICE and anything else required to run them at the time of writing. It also builds a bunch of case directories (under cases) varying in their volumetric flow rate (from 10 m3 to 100 m3). To start a coupled simulation enter the solid and fluid directories of a case with two seperate terminals and execute the appropriate run script.

**Problems during installation:**  
This installation script was built on Ubuntu 22.04.4 LTS (Jammy Jellyfish) running on WSL on Windows 11. Verify the OpenFOAM repo (and others) are compatible or commence an edit in the installPreciceAndOF.sh file. If the installation fails, someone probably updated their version without telling anyone else, making it necessary to fall back to fixed versions compatible with one another.

Everything is intended to run on linux, but this file can also be compiled on windows; just make sure the LF end-of-line flag is set or else all created files will be bogus.

**Just a heads-up:**  
This file contains a string mailPush, which is originally set to my e-mail address. Please keep this in mind before using the generated files to queue a HPC job - or we'll be in touch. Default walltime is set to 4h but will need to be adjusted depending on the available cores and hardware.

**Arguments:**  
imfdFlap takes either "prepare", "build" or "wcgw", which is my version of help because I think it's funny, as an argument. Chances are I ran into your problem setting up precice and openfoam in which case I will explain it there. "imfdFlap build" is automatically invoked by allrun.sh.

**Instructions:**  
Run "./a.out prepare + <name of directory to contain all files (defaults to: imfdFlap)>" to generate installation scripts and start the installation with "allrun.sh".

#### Copy of the "What Could Go Wrong?" section:

Welcome to the What Could Go Wrong? section of my file. If you turn up here, chances are you're miserably stuck so here's a list of issues I encountered that may help you unstick yourself.

If I haven't described the issue you're facing below, I may have wrongly passed it off as a neglible nuisance. Should this be the case, please contact me erstwhile at alrik-matti.luda@student.tu-freiberg.de until I can be bothered to finish setting up github.

**Inaccurate results:**  
Having adapted everything from the 2D perpendicular-flap tutorial found on the preCICE website (source: https://precice.org/tutorials-perpendicular-flap.html), I wasn't suprised to find almost every setting in the solid and fluid fvScheme set to some first order scheme. The only significant improvements in terms of reproducibility (compared to other solvers coupled to openfoam) were achieved by setting

gradSchemes default -> fourth

divSchemes default -> Gauss linear

**Serial run works but parallel run fails with processorPolyPatch error:**  
Take a look at the solution I described here:  https://www.cfd-online.com/Forums/openfoam-solving/255295-pimplefoam-aborts-during-parallel-run-processorpolypatch-error.html

**Serial run works but parallel run results in bad communication:**  
Check the precice-config and find out if you have a loopback device addressed by "lo" or if it has a different address.

**Managing a high Courant number:**  
The Courant number evolves over time as the mesh is moved and reflects the highest velocity in the smallest length scale. Slamming the fluid into the flap without a ramp may squish the cells beyond the point of no return. The solver will abort as soon as the cells crash into each other. Even if the Courant number declines after the initial bounce, you might now be facing a number significantly below 1 with no way to use adjustableTimeStep due to a lack of support for adjustableTimeStep by the openfoam-adapter. Tl;dr: use a ramp (default setting in the generated U BC).

**Coupling is successful but solidDisplacementFoam returns 0 on all calculations:**  
This may be caused by the custom force BC for solidDisplacementFoam not decomposing correctly. I implemented a workaround in the runSolid.sh script so I hope you will never encounter it. Essentially: copy the flap BC from 0/D to all proc\*/0/D's.

**Can't resume simulations from latest time:**  
precice-config.xml needs to be set to a later time as well and the start time adjusted.

**How can I track the progress of my coupled simulation:**  
Besides tailing the solver logs you can take a look at the convergence and iteration logs produced by precice or run the plot-displacement.sh script to visualize your progress so far.

**I want to run the turek-hron-fsi3 case but I don't have groovyBC installed:**  
Glad you asked, here's an installation guide, which wasn't all that easy to find:  
Download and install swak4foam, which is required by the turek-hron-fsi3 precice tutorial  
cd "$HOME/OpenFOAM/$USER-$WM_PROJECT_VERSION"  
hg clone http://hg.code.sf.net/p/openfoam-extend/swak4Foam swak4Foam  
cd swak4Foam  
hg update develop  
echo "Run ./AllwmakeAll -j -q -l inside the openfoam shell"  
openfoam2306  

**I want to try coupling deal.II:**  
Sure, here you go:  
sudo apt install libdeal.ii-dev libdeal.ii-doc cmake make g++  
\#this is optional, build tutorials to test deal.ii  
cp -r /usr/share/doc/libdeal.ii-doc/examples/step-1 .  
cd step-1  
cmake .  
make run  
\#the precice-dealii-adapter  
git clone https://github.com/precice/dealii-adapter.git && cd dealii-adapter  
\#defaults to -DDIM=2 and builds a 2D solver called elasticity  
cmake . //cmake -DDIM=3 .  
make  
\#add the elasticity executables path to your environment or just copy it into the case directory  

**I want to try coupling FEniCs:**  
Sure, here you go:  
sudo apt install software-properties-common  
sudo add-apt-repository ppa:fenics-packages/fenics  
sudo apt update  
sudo apt install fenicsx  
sudo apt install python3-pip  
pip3 install --user fenicsprecice  
\#change lines in perpendicular-flap tutorial case because ufl has been removed from fenics-packages/fenics  
change in the solid.py:  
replace: from ufl import nabla_div  
with: from ufl_legacy import nabla_div  
        import ufl_legacy as ufl  

**I want to try coupling Calculix:**  
Sure, but it will only work with preCICE v2 here you go:  
sudo apt install libarpack2-dev libspooles-dev libyaml-cpp-dev  
wget http://www.dhondt.de/ccx_2.20.src.tar.bz2  
tar xvjf ccx_2.20.src.tar.bz2  
\#make sure you have the right package for your OS, this one is for jammy  
wget https://github.com/precice/calculix-adapter/releases/download/v2.20.0/calculix-precice2_2.20.0-1_amd64_jammy.deb  
sudo apt install ./calculix-precice2_2.20.0-1_amd64_jammy.deb  
cd CalculiX  
wget https://github.com/precice/calculix-adapter/archive/refs/heads/master.tar.gz  
tar -xzf master.tar.gz  
\#fork to take care of a broken file  
cd ccx_2.20/src  
rm cubtri.f  
\#retrieve replacement from https://calculix.discourse.group/t/compilation-error-with-gfortran10/1048/3  
cp /location/of/replacement/cubtri.f .  
cd ../..  
\#end of fork  
cd calculix-adapter-master  
make  
export PATH="/home/ocelittle/CalculiX/calculix-adapter-master/bin:${PATH}"  

**The simulations occupy a lot of disk space:**  
writeCompression is off by default but I string to manage it is available in the main function.  

Where did you find values and equations for your turbulence values:  
Right here and here:  
https://www.openfoam.com/documentation/guides/latest/doc/guide-turbulence-ras-k-omega-sst.html

https://www.cfd-online.com/Wiki/Turbulence_intensity  

These are the most major hurdles I had to take, but there are a lot more to trip over; unfortunately, I forgot to write them down. Add them yourself or message me.



Happy coupling!

-- Alrik Luda, Freiberg, 2024-04-10
