TITLE 'H10.2 ainabea Flex'     { the problem identification }
COORDINATES cartesian3  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
u
v
w
SELECT         { method controls }
!spectral_colors
!ngrid=15
stages = 10

DEFINITIONS    { parameter definitions }
mag = .3*globalmax(magnitude(x,y,z))/globalmax(magnitude(u,v,w))

Lx = .7
Ly = .2
Lz = 3

omega=572.1+.1*stage !determined through trial error

E = if(y<Ly/2) then 100e9 else 50e9
nu = if(y<Ly/2) then 0.3 else 0.5
rho = if(y<Ly/2) then 1000 else 3000
G = E/(2*(1+nu))

ex = dx(u)
ey = dy(v)
ez = dz(w)

gyz = dz(v) + dy(w)
gxz = dz(u) + dx(w)
gxy = dy(u) + dx(v)

C11 = E/((1+nu)*(1-2*nu))*(1-nu)
C12 = E/((1+nu)*(1-2*nu))*nu
C13 =C12
C21 =C12
C22 = C11
C23 = C12
C31 = C12
C32 = C12
C33 = C11

sx = C11*ex + C12*ey+C13*ez
sy = C21*ex + C22*ey+C23*ez
sz = C31*ex + C32*ey+C33*ez

syz = G*gyz
sxz = G*gxz
sxy = G*gxy

phi = atan2(y,x)
xp = x+u
yp = y+v
thetatest = atan2(yp,xp)-phi
theta = if(thetatest<-pi) then thetatest+2*pi else if(thetatest > +pi) then thetatest-2*pi else thetatest

! INITIAL VALUES
EQUATIONS        { PDE's, one for each variable }
u: dx(sx) + dy(sxy) + dz(sxz) = -rho*omega^2*u
v: dx(sxy) + dy(sy) + dz(syz) =  -rho*omega^2*v
w: dx(sxz) + dy(syz) + dz(sz) =  -rho*omega^2*w

EXTRUSION
surface 'bottom' z = 0
surface 'top' z = Lz

BOUNDARIES       { The domain definition }
surface 'bottom'
value(u) = 0
value(v) = 0
load(w) = 0
surface 'top'
value(u) = 0
value(v) = 0
value(w) = 0

!load(w) = 735294.1176*y+73529.41176


  REGION 1       { For each material region }
    START (0,0) line to (Lx,0)
line to (Lx,Ly) 
load(v) = 5  line to (0,Ly)
load(v) = 0	 line to close
	!start(Lx,0) arc(center=0,0) angle=360

! TIME 0 TO 1    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
	grid(x+u*mag, y+v*mag,z+w*mag)
history(w) at (Lx/2,Ly/2, Lz/2) report(omega)
history(v) at (Lx/2,Ly/2, Lz/2) report(omega)
history(u) at (Lx/2,Ly/2, Lz/2) report(omega)
elevation(v) from(Lx/2,Ly/2,0) to (Lx/2,Ly/2,Lz)

SUMMARY
report val(w, Lx,Ly,Lz)

END
