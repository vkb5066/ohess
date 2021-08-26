#!/bin/csh

set origLoc=`pwd`
set dirControl="$origLoc/"
set maxDepth=10
set i=0
set add="*/"

#make dirControl in folder to avoid early termination
set z=0
cd $dirControl
mkdir depthControl
cd depthControl
while ($z < $maxDepth)
        mkdir $z
        cd $z
        @ z++
end


echo "dir,workDir,nrg,nConv,irredKpts,aV,bV,cV,alpha,beta,gamma,vol,sxx,syy,szz,sxy,syz,szx" >$origLoc/output.csv


#Begin looping thorugh sub directories
while ($i < $maxDepth)

	foreach j ($dirControl)
		cd j
	
		if (((-f POSCAR) || (-f CONTCAR)) && -f OUTCAR) then

#do stuff here------------------------------------------------------------------------------

#Find energy...
grep "free  energy" OUTCAR | tail -1 >tmpEnergy
sed -i 's/^.*= \(.*\) eV.*$/\1/p' tmpEnergy
set nrg=`tail -1 tmpEnergy`
rm tmpEnergy

#Find working directory...
set wd=`echo "$j" | awk -F "/" '{print $(NF-1)}'`

#Get number of converged steps
set nConv=`grep "free  energy" OUTCAR | wc -l`

#Get number of irreducible k-points
set ikp=`grep "irreducible" OUTCAR | awk '{print $2}'`

#There should be a CONTCAR but just in case there isnt, use the POSCAR  
if (-f POSCAR) then
	cp POSCAR file__
endif
if (-f CONTCAR) then
	cp CONTCAR file__
endif

#Get the area of the supercell (without C++...oh boy)
set usf=`cat file__ | sed -n '2p'`
set ax=`cat file__ | sed -n '3p' | awk '{print $1}'`
set ay=`cat file__ | sed -n '3p' | awk '{print $2}'`
set az=`cat file__ | sed -n '3p' | awk '{print $3}'`
set bx=`cat file__ | sed -n '4p' | awk '{print $1}'`
set by=`cat file__ | sed -n '4p' | awk '{print $2}'`
set bz=`cat file__ | sed -n '4p' | awk '{print $3}'`
set cx=`cat file__ | sed -n '5p' | awk '{print $1}'`
set cy=`cat file__ | sed -n '5p' | awk '{print $2}'`
set cz=`cat file__ | sed -n '5p' | awk '{print $3}'`	#no trig functions (as far as I know) in bc, so need to explicitly write out the dot of the cross product
set vol=`echo ""$usf"*"$usf"*"$usf"*("$cx"*("$ay"*"$bz" - "$az"*"$by") + "$cy"*("$az"*"$bx" - "$ax"*"$bz") + "$cz"*("$ax"*"$by" - "$ay"*"$bx"))" | bc -l`

#Get the lattice vector lengths for later
set xV=`echo "sqrt("$ax"*"$ax" + "$ay"*"$ay" + "$az"*"$az")" | bc -l`
set yV=`echo "sqrt("$bx"*"$bx" + "$by"*"$by" + "$bz"*"$bz")" | bc -l`
set zV=`echo "sqrt("$cx"*"$cx" + "$cy"*"$cy" + "$cz"*"$cz")" | bc -l`

#Get cosines of the angles
set cosAl=`echo "(("$bx"*"$cx")+("$by"*"$cy")+("$bz"*"$cz"))/("$yV"*"$zV")" | bc -l`
set cosBe=`echo "(("$ax"*"$cx")+("$ay"*"$cy")+("$az"*"$cz"))/("$xV"*"$zV")" | bc -l`
set cosGa=`echo "(("$bx"*"$ax")+("$by"*"$ay")+("$bz"*"$az"))/("$yV"*"$xV")" | bc -l`

#Get Angles - turns out there are limited trig functions in bc
##alpha
if (`echo ""$cosAl" == 0" | bc -l`) then
       set al="90"
else if (`echo ""$cosAl" <= 1" | bc -l`) then
       set al=`echo "a(sqrt((1/("$cosAl"*"$cosAl"))-1))*180/3.14159" | bc -l`
else
	set al="ERR"
endif
##beta
if ((`echo ""$cosBe" == 0" | bc -l`)) then
       set be="90"
else if ((`echo ""$cosBe" <= 1" | bc -l`)) then
       set be=`echo "a(sqrt((1/("$cosBe"*"$cosBe"))-1))*180/3.14159" | bc -l`
else
	set be="ERR"
endif
##gamma
if ((`echo ""$cosGa" == 0" | bc -l`)) then
       set ga="90"
else if ((`echo ""$cosGa" <= 1" | bc -l`)) then
       set ga=`echo "a(sqrt((1/("$cosGa"*"$cosGa"))-1))*180/3.14159" | bc -l`
else
	set ga="ERR"
endif

#Get Stresses
set strsXX=`grep "in kB" OUTCAR | tail -1 | awk '{print $3}'`
set strsYY=`grep "in kB" OUTCAR | tail -1 | awk '{print $4}'`
set strsZZ=`grep "in kB" OUTCAR | tail -1 | awk '{print $5}'`
set strsXY=`grep "in kB" OUTCAR | tail -1 | awk '{print $6}'`
set strsYZ=`grep "in kB" OUTCAR | tail -1 | awk '{print $7}'`
set strsZX=`grep "in kB" OUTCAR | tail -1 | awk '{print $8}'`

#Append all of this information to the csv line
echo "`pwd`,$wd,$nrg,$nConv,$ikp,$xV,$yV,$zV,$al,$be,$ga,$vol,$strsXX,$strsYY,$strsZZ,$strsXY,$strsYZ,$strsZX" >>$origLoc/output.csv

rm file__

#stop doing stuff here----------------------------------------------------------------------
		
		cd ../
		endif
	end

set dirControl="$dirControl$add"
@ i++

end

cd $origLoc
rm -rf depthControl
echo "Script End"

