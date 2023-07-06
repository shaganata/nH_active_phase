------------------------------------------------------------------------

  INSTRUCTIONS TO USE THE SOFTWARE FOR FITTING H^0 COLUMN DENSITIES
    OF HIGH-INCLINATION SYMBIOTIC STARS DURING THEIR ACTIVE PHASES
        
------------------------------------------------------------------------  
CONTENT:


FOLDER  "/code/":
        
01_nH_inversion_fit_5_2.f90    - code file for egress/ingress column 
                                density fitting and velocity profile 
                                derivation, further referenced as 
                                "code 01"   
                                
02_nH_act_phase_all_5_5.f90    - code file for computation of column 
                                density at all orbital phases, further 
                                referenced as "code 02", combine the 
                                results of fitting by code 01                                
                                       
FOLDER  "/data/":

   SUBFOLDER: "/01_nH_egr_ingr_inversion/":

   all-nh_fi_EGRESS.dat          - test data file, further referenced as 
                                   "data 01a"     
   all-nh_fi_INGRESS.dat         - test data file, further referenced as 
                                   "data 01b"                            
   eigenvaluesAbel.dat           - eigenvalues of Abel operator to be read 
                                   by code 01
   fit_Egress_i90_nHwd1p5e22.txt - expected output of code 01 after several
                                   runs, always with narrowed ranges of 
                                   possible values of fitting parameters                            
        
   SUBFOLDER  "/02_nH_active_phase_all/":
   
   nH_act_phase_measured.txt     - test data file, further referenced as 
                                   "data 02"
   nH_act_phase_i90_wd1c5e22.txt - expected output of the code 02                            
   nH_act_phase.gnu              - gnuplot script for data and model 
                                   visualisation
   nH_act_phase_i90_wd1c5e22.eps - expected output of the gnuplot script

------------------------------------------------------------------------
1. SYSTEM REQUIREMENTS:

Hardware requirements: standard computer

OS requirements: Linux, Windows or MacOS
                The package has been tested on the following systems:
                Linux Mint 18.3, Centos 7, Windows 10 Pro v1909

Software requirements: Fortran 90 compiler (for the list of compilers, 
                       see: https://fortran-lang.org/compilers/).
                       
                       Gnuplot (for data visualisation with
                       nH_act_phase.gnu script) or other data plotting
                       software of your choise.

------------------------------------------------------------------------
2. INSTALLATION GUIDE:

Fortran 90 compiler is required. For example, instructions to install 
GFortran can be found at:
https://fortran-lang.org/learn/os_setup/install_gfortran

For Windows systems, GFortran can be used from within the Msys2 platform,
see the installation instructions at:
https://masuday.github.io/fortran_tutorial/install_gfortran_windows.html

No other installation for codes 01 and 02 is required.

Typical install time: ~ minutes
------------------------------------------------------------------------
3. DEMO:

INSTRUCTIONS TO RUN ON DATA

The codes 01 and 02 are written in Fortran 90 language and can be 
launched in the command line by commands depending on the Fortran 
compiler used and the operating system.

Examples of commands to compile and run the codes from withing the 
directory with /code/ and /data/ folders:

....................

Linux system and GFortran compiler:

cd ./code
gfortran 01_nH_inversion_fit_5_2.f90 -o 01_nH_inversion_fit_5_2.out
cd ../data/01_nH_egr_ingr_inversion/
../../code/01_nH_inversion_fit_5_2.out

cd ./code
gfortran 02_nH_act_phase_all_5_5.f90 -o 02_nH_act_phase_all_5_5.out
cd ../data/02_nH_active_phase_all/
../../code/02_nH_act_phase_all_5_5.out

....................

Windows system and GFortran compiler using the bash command line 
terminal within Msys2 platform: 

cd ./code
gfortran 01_nH_inversion_fit_5_2.f90 -o 01_nH_inversion_fit_5_2.exe
cd ../data/01_nH_egr_ingr_inversion/
../../code/01_nH_inversion_fit_5_2.exe

cd ./code
gfortran 02_nH_act_phase_all_5_5.f90 -o 02_nH_act_phase_all_5_5.exe
cd ../data/02_nH_active_phase_all/
../../code/02_nH_act_phase_all_5_5.exe

.......................................................................

Running time on a standard computer:

code 01:  1 - 30 minutes for a dataset comparable with data 01a or 01b 
          and 50000 iterations 
code 02:  ~ few seconds

.......................................................................

The gnuplot script nH_act_phase.gnu can be run in the command line by 
command:

gnuplot nH_act_phase.gnu

The output .eps file is then created.

------------------------------------------------------------------------
4. INSTRUCTIONS FOR USE: 

MODEL:

....................

code 01:

Egress/ingress column density originating in the giant wind is fitted by
a polynomial expression: 
n_H(b) = n1*b^{-1} + nK*b^{-K}, 
where n1, nK and K are parameters, b is impact parameter
(Dumm et al. 1999, A&A 349, 169).
The contribution from the hot component, nHwd, is assumed to be 
constant.

The definition of impact parameter b allows to fit the data observed at 
the orbital phases <0,0.25> (egress) or <0.75,1> (ingress). At this 
orbital phases, the nHwd value is low in comparison with n_H(b) coming
from giant wind. Therefore, its value is determined by using the code 01
together with the code 02.

Corresponding velocity profile is obtained by inversion of an Abel-type 
operator for column density function (Knill et al. 1993, A&A 274, 1002)
in the form:
v(r)=v_ter/(1+xi*r^{1-K})
where r is the radial distance from the giant centre, v_ter is terminal 
velocity, xi=nK*lambda1/(n1*lambdaK), lambda1 and lambdaK are 
eigenvalues of Abel operator.

The value of n1 parameter determines the Mdot/v_ter value, where Mdot is
mass-loss rate from the giant (Shagatova et al. 2016, A&A 588, A83).

....................

code 02:

The model enables interconnection of two velocity profiles in the case 
of asymmetric density distributions in the plane of observation through
"transitional" velocity profile, assuming smooth change of the velocity 
profile from the egress one to the ingress profile in a region defined 
by a boundary values of the angle psi (Shagatova et al. 2017, A&A 602, 
A71). In this way, the column density of the giant wind is computed at 
all orbital phases using the continuity equation.

This code enables optimalization of the nHwd value, which contributes 
significantly to the total measured column density at phi ~ 0.25 - 0.75.
Reduced Chi^2 value is given for the complete dataset from egress to 
ingress orbital phases.

.......................................................................

FITTING PROCEDURE:

Part I.: Fitting the egress/ingress observed column densities (i.e. 
orbital phase phi <=0.25 / >=0.75):
1. Enter the input to the code 01: data to be fitted (e.g. data 01a, 
   line 59), orbital inclination i (line 60) and the fixed value of 
   column density contribution from the hot component nHwd that will be 
   substracted from the observed neutral hydrogen column density values 
   (line 61).
2. If needed, change the number of iterations to be done (line 159) and 
   the ranges of fitting parameters (lines 281 - 286).
3. Run the code 01 and check the output file with 50 best fits found. 
   If the precision of fitting parameter values is not satisfactory,
   repeat the steps 2.-3. until the acceptable precision is achieved.
4. If you wish to fit another dataset, e.g. from the complementary 
   fraction of orbital phases (ingress/egress), go back to the steps 
   1. - 3.
Result: egress/ingress column density and velocity profile 
   
Part II.: Computation of column density of giant wind at all orbital 
phases:
1. Enter the input to the code 02: 
   - fitting parameters of egress and ingress column density and 
     velocity profile obtained in Part I. (lines 60 - 67)
   - corresponding orbital inclination and nHwd value (lines 57 and 71)
   - the name of output file (line 73)
   - you can change the extent of the transitional region between areas 
     of egress and ingress velocity profile by changing the boundary 
     values of angle psi: psiE and psiIN (lines 79 - 80)
2. Run the code 02 and verify, whether the value of nHwd selected in the 
   step 1. of Part I. match the data (e.g. data 02) at orbital phases 
   0.25 - 0.75 (the value of the reduced Chi^2 is shown in the command 
   line, you can use also gnuplot script nH_act_phase.gnu for 
   visualisation). If not, go back to the step 1. of Part I. and start 
   with a different value of nHwd. Repeat the whole procedure until the 
   adequate value of nHwd is found.
Result: column density as a function of phi originating from the giant
        wind + the value of constant contribution from the hot component
   
.......................................................................

INSTRUCTIONS TO REPRODUCE THE RESULTS IN THE MANUSCRIPT:

Models in the manuscript NCOMMS-20-47961-T were obtained for three 
different orbital inclinations: 70°, 80° and 90°. To reproduce this 
models, the following steps have to be carried out for each value of 
orbital inclination:

  i. Perform steps 1.-3. of Part I. for data 01a to obtain egress column 
     density and velocity profile.
 ii. Perform steps 1.-3. of Part I. for data 01b to obtain ingress 
     column density and velocity profile. 
     Due to the missing ingress data points around orbital phase ~ 0.8, 
     the fitting procedure does not yield the adequate value of n1 
     parameter for orbital inclinations 70° and 80° for the 
     representation of the whole dataset (data 02, see step iii. below). 
     In these cases, the appropriate n1 of ingress model was obtained by 
     fixing its value at 2.2e+23 and performing steps 1.-3. of Part I. 
     for inclination of 70°; and fixing n1=2.05e+23* and K=5 for 
     inclination of 80°. For inclination of 90°, no manual intervention 
     is needed.
iii. Perform steps 1.-2. of Part II. for data 02 to obtain the global 
     column density model. The models in the manuscript correspond to 
     the nHwd values 0.5e+22, 1.0e+22 and 1.5e+22 for inclinations of 
     70°, 80° and 90°, respectively. Due to the high uncertainty of this
     value, one can arrive to different nHwd values by repeating steps 
     in Part I. and II., but generally values around 1.0e+22 provide an
     acceptable model.
     
* the ingress n1(80°) value was estimated from the ratio of obtained 
  egress n1 values:
  (n1(70°)-n1(80°))/(n1(80°)-n1(90°)),
  by applying the same factor to the ingress n1 values.
  
------------------------------------------------------------------------
5. SOFTWARE LICENCE:

GNU General Public License v3.0

https://opensource.org/licenses/GPL-3.0

------------------------------------------------------------------------
