from yade import ymport, utils, plot

#################################### parameters definition
DAMPING=0.5
dtCOEFF=0.5
normalSTRESS=150e3
normalVEL=normalSTRESS/1e7 
shearVEL=1*normalVEL 
intR=1.263 
DENS=2700
YOUNG1=150e6
YOUNG2=50e9
FRICT=18
ALPHA=0.3
TENS=45e5
COH=45e6
iterMAX=400000

###################### Import of the sphere assembly
wallMat=O.materials.append(FrictMat(young=YOUNG2,frictionAngle=0,density=DENS,poisson=.1,label="wallMat"))

aggregate=O.materials.append(FrictMat(density=2.6e3, young=150e6,poisson=.2,frictionAngle=0.785,label="aggregate"))
ymport.textClumps('compress.txt',scale=1.,shift=Vector3(0,0,0),material=aggregate)
N1=O.bodies[-1].id

## preprocessing to get dimensions
dim=utils.aabbExtrema()
xinf=dim[0][0]
xsup=dim[1][0]
X=xsup-xinf
yinf=dim[0][1]
ysup=dim[1][1]
Y=ysup-yinf
zinf=dim[0][2]
zsup=dim[1][2]
Z=zsup-zinf
H=0.147

S0=X*Y


### loading platens
O.bodies.append(utils.box(center=(xinf+X/2.,yinf+Y/2.,zinf-0.005),extents=(X/2,Y/2, 0.005),material=wallMat,fixed=True))
bottomPlate=O.bodies[-1]
O.bodies.append(utils.box(center=(xinf+X/2.,yinf+Y/2.,zsup+0.005),extents=(X/2,Y/2, 0.005),material=wallMat,fixed=True))
topPlate=O.bodies[-1]


O.bodies.append(utils.box(center=(xinf-X/4-0.01, yinf+Y/2.,H-0.005),extents=(X/4,Y/2, 0.005),material=wallMat,fixed=True))
WALL=O.bodies[-1]

### add top box
O.bodies.append(utils.box(center=(xinf-0.005,0,0.26),extents=(0.005,Y/2+0.02,0.1),material=wallMat,color=(0.5, 0.5, 0.5), fixed=True,wire=True))
face1=O.bodies[-1]
O.bodies.append(utils.box(center=(xsup+0.005,0,0.26),extents=(0.005,Y/2+0.02,0.1),material=wallMat,color=(0.5, 0.5, 0.5), fixed=True,wire=True))
face2=O.bodies[-1]
O.bodies.append(utils.box(center=(0, yinf-0.005,0.26),extents=(X/2,0.005,0.1),material=wallMat,color=(0.5, 0.5, 0.5),fixed=True,wire=True))
face3=O.bodies[-1]
O.bodies.append(utils.box(center=(0, ysup+0.005,0.26),extents=(X/2,0.005,0.1),material=wallMat,color=(0.5, 0.5, 0.5),fixed=True,wire=True))
face4=O.bodies[-1]


### add bottom box
O.bodies.append(utils.box(center=(xinf-0.005,0,H/2),extents=(0.005,Y/2+0.02,H/2),material=wallMat,color=(0.5, 0.5, 0.5), fixed=True,wire=True))
face5=O.bodies[-1]
O.bodies.append(utils.box(center=(xsup+0.005,0,H/2),extents=(0.005,Y/2+0.02,H/2),material=wallMat,color=(0.5, 0.5, 0.5), fixed=True,wire=True))
face6=O.bodies[-1]
O.bodies.append(utils.box(center=(0, yinf-0.005,H/2),extents=(X/2,0.005,H/2),material=wallMat,color=(0.5, 0.5, 0.5),fixed=True,wire=True))
face7=O.bodies[-1]
O.bodies.append(utils.box(center=(0, ysup+0.005,H/2),extents=(X/2,0.005,H/2),material=wallMat,color=(0.5, 0.5, 0.5),fixed=True,wire=True))
face8=O.bodies[-1]



O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb(),Bo1_Box_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom_FrictPhys_CundallStrack()]
	),
	VTKRecorder(iterPeriod=200000,recorders=['intr'],fileName='chain'),
	TranslationEngine(ids=[topPlate.id],translationAxis=(0.,0,-1),velocity=0.,label='zTranslation'),
	NewtonIntegrator(damping=DAMPING,gravity=(0,0,-9.8),label='damper'),
	PyRunner(command='shearLoading()',iterPeriod=3,label='checker'),
	PyRunner(command='stopShearing()',iterPeriod=1000,label='checker'),
	PyRunner(command='addPlotData()',iterPeriod=2000),
]

O.dt=2e-6


shearing=False
tau=0
Fs1=0
Fs2=0
Xdispl=0
px0=0
Zdispl=0
pz0=topPlate.state.pos[2]
prevTranslation=0
n=0


def shearLoading():
	global px0, pz0, n, shearing, sigmaN
	px0=face6.state.pos[0]-0.01
	pz0=topPlate.state.pos[2]
	sigmaN=O.forces.f(topPlate.id)[2]/S0
	
	if shearing==False:
		if zTranslation.velocity<normalVEL:
			zTranslation.velocity+=normalVEL/1000
		if sigmaN>(0.975*normalSTRESS):
			zTranslation.velocity=normalVEL*((normalSTRESS-sigmaN)/normalSTRESS)

	if shearing==False and abs((normalSTRESS-sigmaN)/normalSTRESS)<0.025:
		zTranslation.velocity=0
		n+=1
		if n>20 and abs((sigmaN-normalSTRESS)/normalSTRESS)<0.025:
			O.engines=O.engines+ [TranslationEngine(ids=[face5.id,face6.id, face7.id, face8.id, bottomPlate.id, WALL.id],translationAxis=(1.,0.,0.),velocity=0.001)]

#			px0=face5.state.pos[0]
#			pz0=topPlate.state.pos[2]
			shearing=True
			print ('shearing now! || iteration=', O.iter)
			O.save('0.yade.gz')
	if shearing==True:
		zTranslation.velocity=200*normalVEL*((normalSTRESS-O.forces.f(topPlate.id)[2]/S0)/normalSTRESS)
                 

def stopShearing():
	if abs(face6.state.pos[0]-X/2-0.005)/X > 0.005 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.00502:
		O.save('1.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.01 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.01002:
		O.save('2.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.015 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.01502:
		O.save('3.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.02 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.02002:
		O.save('4.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.025 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.02502:
		O.save('5.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.03 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.03002:
		O.save('6.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.035 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.03502:
		O.save('7.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.04 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.04002:
		O.save('8.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.045 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.04502:
		O.save('9.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.05 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.05002:
		O.save('10.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.055 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.05502:
		O.save('11.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.06 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.06002:
		O.save('12.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.065 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.06502:
		O.save('13.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.07 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.07002:
		O.save('14.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.075 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.07502:
		O.save('15.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.08 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.08002:
		O.save('16.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.085 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.08502:
		O.save('17.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.09 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.09002:
		O.save('18.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.095 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.09502:
		O.save('19.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.10 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.10002:
		O.save('20.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.105 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.10502:
		O.save('21.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.11 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.11002:
		O.save('22.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.115 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.11502:
		O.save('23.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.12 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.12002:
		O.save('24.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.125 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.12502:
		O.save('25.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.13 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.13002:
		O.save('26.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.135 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.13502:
		O.save('27.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.14 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.14002:
		O.save('28.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.145 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.14502:
		O.save('29.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.15 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.15002:
		O.save('30.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.155 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.15502:
		O.save('31.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.16 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.16002:
		O.save('32.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.165 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.16502:
		O.save('33.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.17 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.17002:
		O.save('34.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.175 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.17502:
		O.save('35.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.18 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.18002:
		O.save('36.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.185 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.18502:
		O.save('37.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.19 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.19002:
		O.save('38.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.195 and abs(face6.state.pos[0]-X/2-0.005)/X < 0.19502:
		O.save('39.yade.gz')
	elif abs(face6.state.pos[0]-X/2-0.005)/X > 0.19999:
		plot.saveDataTxt('2.txt')
		O.pause()
		O.save('40.yade.gz')  

             
def addPlotData():
	w=face5.state.pos[0]
	Fx=O.forces.f(face5.id)[0]
	Fz=O.forces.f(topPlate.id)[2]/(0.5*(0.5-abs(face6.state.pos[0]-X/2-0.005)))
	zdispl=topPlate.state.pos[2]
	CN=avgNumInteractions(cutoff=0.,skipFree=False,considerClumps=True)
	plot.addData(i=O.iter, j=O.iter, k=O.iter, Fx=Fx, Fz=Fz, w=w,zdispl=zdispl, CN=CN)


plot.plots={'w': ('Fx'),'i': ('zdispl'), 'j': ('Fz'), 'k': ('CN')}
plot.plot()

