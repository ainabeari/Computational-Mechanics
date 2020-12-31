TITLE 'H8.2 ainabea FlexPDE'
COORDINATES cartesian3
VARIABLES        { system variables }
	u   
	v
	w
SELECT         { method controls }
ngrid=20
DEFINITIONS    { parameter definitions }
mag=1e3

grav=9.81

Lz = 10

!piecewise functions
E=  if x<0 then 210e9 else 5e9
nu= if x<0 then 0.4 else 0.2
rho= if x<0 then 7600 else 990

G=E/(2*(1+nu))

C11 =E*(1-nu)/((1+nu)*(1-2*nu))
C22 = C11
C33 = C11

C12 = E*nu/((1+nu)*(1-2*nu))
C13 = C12
C21 = C12
C23 = C12
C31 = C12
C32 = C12

!!Strain
!Axial Strain
ex=dx(u)
ey=dy(v)
ez=dz(w)
!Engineering Shear Strain
gxy=dx(v)+dy(u)
gyz=dy(w)+dz(v)
gxz=dz(u)+dx(w)

!Stress via Hooke's Laws
!Axial Stress
sx = C11*ex + C12*ey + C13*ez
sy = C21*ex + C22*ey + C23*ez
sz = C31*ex + C32*ey + C33*ez
!Shear Stress
syz= G*gyz
sxz= G*gxz 
sxy= G*gxy

EQUATIONS        { PDE's, one for each variable }
!Fnet = 0
u:		dx(sx)+dy(sxy)+dz(sxz)=0
v:		dx(sxy)+dy(sy)+dz(syz)-rho*grav + z^2*10000=0
w:	dx(sxz)+dy(syz)+dz(sz)=0

EXTRUSION
surface 'bottom' z=0
surface 'top' z=Lz

BOUNDARIES       { The domain definition }
surface 'bottom'
	value(u)=0
	value(v)=0
	value(w)=0
surface 'top' !remember to say boundary conditions!
	!value(u)=0
	!value(v)=0
	!value(w)=0

  REGION 1       { For each material region }
    START(-0.25,0)
	load(u)=0
	load(v)=0
	load(w)=0
	LINE TO(0.25,0)TO(0.25,0.2)TO(0.05,0.2)TO(0.05,0.6)TO(0.25,0.6)TO(0.25,0.8)TO(-0.25,0.8)TO(-0.25,0.6)TO(-0.05,0.6)TO(-0.05,0.2)TO(-0.25,0.2)TO CLOSE
PLOTS            
	grid(x+mag*u,y+mag*v,z+mag*w)
	contour(sz) painted on z=Lz/2
	contour(ez) painted on z=Lz/2
	elevation(v) from (0,0.55,0) to (0,0.55,Lz)
SUMMARY
	report val(v,0,0,Lz/2)
	report val(v,0,0.55,Lz/2)
END

