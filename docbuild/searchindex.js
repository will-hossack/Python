Search.setIndex({docnames:["analysis","index","jones","lens","material","matrix","psf","ray","surface","tio","vector","wavefront","wavelength","zernike"],envversion:{"sphinx.domains.c":1,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":1,"sphinx.domains.index":1,"sphinx.domains.javascript":1,"sphinx.domains.math":2,"sphinx.domains.python":1,"sphinx.domains.rst":1,"sphinx.domains.std":1,"sphinx.ext.intersphinx":1,"sphinx.ext.todo":2,"sphinx.ext.viewcode":1,sphinx:56},filenames:["analysis.rst","index.rst","jones.rst","lens.rst","material.rst","matrix.rst","psf.rst","ray.rst","surface.rst","tio.rst","vector.rst","wavefront.rst","wavelength.rst","zernike.rst"],objects:{"optics.analysis":{KnifeTest:[0,0,1,""],OpticalImage:[0,0,1,""],TargetPlane:[0,0,1,""]},"optics.analysis.KnifeTest":{getImage:[0,1,1,""],setKnife:[0,1,1,""],setWire:[0,1,1,""]},"optics.analysis.OpticalImage":{addTestGrid:[0,1,1,""],draw:[0,1,1,""],getImage:[0,1,1,""],getPixelSourcePoint:[0,1,1,""],getRayPencil:[0,1,1,""],getSurfaceInteraction:[0,1,1,""],getSystemImage:[0,1,1,""]},"optics.analysis.TargetPlane":{add:[0,1,1,""],addGrid:[0,1,1,""],draw:[0,1,1,""],getPencils:[0,1,1,""],rayPencil:[0,1,1,""]},"optics.jones":{HalfWavePlate:[2,0,1,""],JonesMatrix:[2,0,1,""],JonesMatrixSystem:[2,0,1,""],JonesVector:[2,0,1,""],LeftCircularPolarisedBeam:[2,0,1,""],LinearPolarisedBeam:[2,0,1,""],LinearPolariser:[2,0,1,""],QuarterWavePlate:[2,0,1,""],Retarder:[2,0,1,""],RightCircularPolarisedBeam:[2,0,1,""]},"optics.jones.JonesMatrix":{copy:[2,1,1,""],determinant:[2,1,1,""],rotate:[2,1,1,""],rotateBy:[2,1,1,""],set:[2,1,1,""],setAngle:[2,1,1,""],setWavelength:[2,1,1,""],trace:[2,1,1,""]},"optics.jones.JonesMatrixSystem":{getMatrix:[2,1,1,""],polarPlot:[2,1,1,""],rotateBy:[2,1,1,""],setWavelength:[2,1,1,""]},"optics.jones.JonesVector":{copy:[2,1,1,""],getAngle:[2,1,1,""],getEllipicity:[2,1,1,""],getIntensity:[2,1,1,""],getPhase:[2,1,1,""],polarPlot:[2,1,1,""],throughPolariser:[2,1,1,""]},"optics.jones.LinearPolariser":{copy:[2,1,1,""],rotate:[2,1,1,""],rotateBy:[2,1,1,""],setAngle:[2,1,1,""]},"optics.jones.Retarder":{copy:[2,1,1,""],setAngle:[2,1,1,""]},"optics.lens":{DataBaseLens:[3,0,1,""],Doublet:[3,0,1,""],Lens:[3,0,1,""],OpticalGroup:[3,0,1,""],SimpleSinglet:[3,0,1,""],Singlet:[3,0,1,""],getCurrentLens:[3,2,1,""],setCurrentLens:[3,2,1,""]},"optics.lens.Doublet":{draw:[3,1,1,""],getFNo:[3,1,1,""],getRadius:[3,1,1,""],getSingletPair:[3,1,1,""],invert:[3,1,1,""],setCurvatures:[3,1,1,""],setFocalLength:[3,1,1,""],setRadius:[3,1,1,""]},"optics.lens.Lens":{backFocalLength:[3,1,1,""],backFocalPlane:[3,1,1,""],backNodalPoint:[3,1,1,""],backPrincipalPlane:[3,1,1,""],cardinalPoints:[3,1,1,""],draw:[3,1,1,""],entrancePupil:[3,1,1,""],exitPupil:[3,1,1,""],frontFocalLength:[3,1,1,""],frontFocalPlane:[3,1,1,""],frontNodalPoint:[3,1,1,""],frontPrincipalPlane:[3,1,1,""],petzvalSum:[3,1,1,""],setFocalLength:[3,1,1,""],setIris:[3,1,1,""]},"optics.lens.OpticalGroup":{add:[3,1,1,""],draw:[3,1,1,""],entranceAperture:[3,1,1,""],exitAperture:[3,1,1,""],getInfo:[3,1,1,""],getPoint:[3,1,1,""],imagePoint:[3,1,1,""],movePoint:[3,1,1,""],paraxialGroup:[3,1,1,""],paraxialMatrix:[3,1,1,""],planePair:[3,1,1,""],scale:[3,1,1,""],setPoint:[3,1,1,""],setWithPlanes:[3,1,1,""]},"optics.lens.Singlet":{draw:[3,1,1,""],getBend:[3,1,1,""],getEdgeThickness:[3,1,1,""],getFNo:[3,1,1,""],getRadius:[3,1,1,""],getThickness:[3,1,1,""],invert:[3,1,1,""],setBend:[3,1,1,""],setCentreThickness:[3,1,1,""],setCurvatures:[3,1,1,""],setEdgeThickness:[3,1,1,""],setFocalLength:[3,1,1,""],setFromString:[3,1,1,""],setParameters:[3,1,1,""],setRadius:[3,1,1,""],setThickness:[3,1,1,""]},"optics.material":{Material:[4,0,1,""],MaterialData:[4,0,1,""]},"optics.material.MaterialData":{getList:[4,1,1,""],getMaterial:[4,1,1,""]},"optics.matrix":{CavityMatrix:[5,0,1,""],DataBaseMatrix:[5,0,1,""],DielectricMatrix:[5,0,1,""],DoubletMatrix:[5,0,1,""],MirrorMatrix:[5,0,1,""],ParaxialDoublet:[5,0,1,""],ParaxialGroup:[5,0,1,""],ParaxialMatrix:[5,0,1,""],ParaxialMirror:[5,0,1,""],ParaxialPlane:[5,0,1,""],ParaxialThickLens:[5,0,1,""],ParaxialThinLens:[5,0,1,""],PropagationMatrix:[5,0,1,""],ThickLensMatrix:[5,0,1,""],ThinLensMatrix:[5,0,1,""]},"optics.matrix.ParaxialGroup":{backFocalPlane:[5,1,1,""],backNodalPoint:[5,1,1,""],backPrincipalPlane:[5,1,1,""],cardinalPoints:[5,1,1,""],copy:[5,1,1,""],draw:[5,1,1,""],frontFocalPlane:[5,1,1,""],frontNodalPoint:[5,1,1,""],frontPrincipalPlane:[5,1,1,""],getInfo:[5,1,1,""],getPoint:[5,1,1,""],imagePlane:[5,1,1,""],imagePoint:[5,1,1,""],incrementInputPlane:[5,1,1,""],inputPlane:[5,1,1,""],maxRadius:[5,1,1,""],outputPlane:[5,1,1,""],planePair:[5,1,1,""],scale:[5,1,1,""],setInputPlane:[5,1,1,""],setWithPlanes:[5,1,1,""]},"optics.matrix.ParaxialMatrix":{backFocalLength:[5,1,1,""],backFocalPlane:[5,1,1,""],backPower:[5,1,1,""],backPrincipalPlane:[5,1,1,""],copy:[5,1,1,""],determinant:[5,1,1,""],frontFocalLength:[5,1,1,""],frontFocalPlane:[5,1,1,""],frontPower:[5,1,1,""],frontPrincipalPlane:[5,1,1,""],inverse:[5,1,1,""],scale:[5,1,1,""],setFocalLength:[5,1,1,""],trace:[5,1,1,""]},"optics.matrix.ParaxialPlane":{draw:[5,1,1,""],getInfo:[5,1,1,""]},"optics.psf":{FixedMoments:[6,0,1,""],Psf:[6,0,1,""],SpotDiagram:[6,0,1,""]},"optics.psf.FixedMoments":{addPoint:[6,1,1,""],area:[6,1,1,""],centroid:[6,1,1,""],eccentricity:[6,1,1,""],ellipse:[6,1,1,""],radius:[6,1,1,""]},"optics.psf.Psf":{area:[6,1,1,""],draw:[6,1,1,""],eccentricity:[6,1,1,""],ellipse:[6,1,1,""],optimalArea:[6,1,1,""],setWithRays:[6,1,1,""]},"optics.psf.SpotDiagram":{draw:[6,1,1,""]},"optics.ray":{IntensityRay:[7,0,1,""],ParaxialRay:[7,0,1,""],PrintPath:[7,0,1,""],Ray:[7,0,1,""],RayMonitor:[7,0,1,""],RayPath:[7,0,1,""],RayPencil:[7,0,1,""],SourcePoint:[7,0,1,""]},"optics.ray.IntensityRay":{copy:[7,1,1,""],getPhaselength:[7,1,1,""],isValid:[7,1,1,""],pointInPlane:[7,1,1,""],propagate:[7,1,1,""],propagateThrough:[7,1,1,""],setInvalid:[7,1,1,""]},"optics.ray.ParaxialRay":{copy:[7,1,1,""],crosses:[7,1,1,""],crossesZero:[7,1,1,""],isValid:[7,1,1,""],mult:[7,1,1,""],multBy:[7,1,1,""],propagate:[7,1,1,""],propagateThrough:[7,1,1,""],propagateTo:[7,1,1,""],setInvalid:[7,1,1,""]},"optics.ray.PrintPath":{update:[7,1,1,""]},"optics.ray.Ray":{addMonitor:[7,1,1,""],draw:[7,1,1,""],isValid:[7,1,1,""],setInvalid:[7,1,1,""],updateMonitor:[7,1,1,""]},"optics.ray.RayMonitor":{draw:[7,1,1,""],update:[7,1,1,""]},"optics.ray.RayPath":{copy:[7,1,1,""],draw:[7,1,1,""],getInfo:[7,1,1,""],update:[7,1,1,""]},"optics.ray.RayPencil":{addBeam:[7,1,1,""],addCollimatedBeam:[7,1,1,""],addCollimatedParaxialBeam:[7,1,1,""],addMonitor:[7,1,1,""],addSourceBeam:[7,1,1,""],addSourceParaxialBeam:[7,1,1,""],draw:[7,1,1,""],propagate:[7,1,1,""],propagateThrough:[7,1,1,""],removeInvalid:[7,1,1,""]},"optics.ray.SourcePoint":{copy:[7,1,1,""],getIntensity:[7,1,1,""]},"optics.surface":{AnnularAperture:[8,0,1,""],CircularAperture:[8,0,1,""],Clear:[8,4,1,""],FlatSurface:[8,0,1,""],ImagePlane:[8,0,1,""],IrisAperture:[8,0,1,""],KnifeAperture:[8,0,1,""],OpticalPlane:[8,0,1,""],ParabolicSurface:[8,0,1,""],QuadricSurface:[8,0,1,""],Reflecting:[8,4,1,""],Refracting:[8,4,1,""],SphericalSurface:[8,0,1,""],Surface:[8,0,1,""]},"optics.surface.AnnularAperture":{draw:[8,1,1,""],getNormal:[8,1,1,""],getSurfaceInteraction:[8,1,1,""],scale:[8,1,1,""]},"optics.surface.CircularAperture":{draw:[8,1,1,""],getNormal:[8,1,1,""],getRadius:[8,1,1,""],getSurfaceInteraction:[8,1,1,""],scale:[8,1,1,""]},"optics.surface.FlatSurface":{getDistance:[8,1,1,""],getNormal:[8,1,1,""],getSurfaceInteraction:[8,1,1,""]},"optics.surface.ImagePlane":{draw:[8,1,1,""],getSize:[8,1,1,""],scale:[8,1,1,""],setSize:[8,1,1,""]},"optics.surface.IrisAperture":{draw:[8,1,1,""],getNormal:[8,1,1,""],getRadius:[8,1,1,""],getSurfaceInteraction:[8,1,1,""],setRatio:[8,1,1,""]},"optics.surface.KnifeAperture":{getSurfaceInteraction:[8,1,1,""],setKnife:[8,1,1,""],setWire:[8,1,1,""]},"optics.surface.OpticalPlane":{draw:[8,1,1,""],getSourcePoint:[8,1,1,""],getSurfaceInteraction:[8,1,1,""],incrementSurface:[8,1,1,""],paraxialGroup:[8,1,1,""],surfaceVector:[8,1,1,""]},"optics.surface.QuadricSurface":{draw:[8,1,1,""],eccentricity:[8,1,1,""],edgePlane:[8,1,1,""],getDistance:[8,1,1,""],getNormal:[8,1,1,""],getParaxialInteraction:[8,1,1,""],getSurfaceInteraction:[8,1,1,""],scale:[8,1,1,""]},"optics.surface.Surface":{draw:[8,1,1,""],getNormal:[8,1,1,""],getPoint:[8,1,1,""],getSurfaceInteraction:[8,1,1,""],makeStandAlone:[8,1,1,""],movePoint:[8,1,1,""],scale:[8,1,1,""],setPoint:[8,1,1,""]},"optics.wavefront":{Interferometer:[11,0,1,""],KingslakeWaveFront:[11,0,1,""],WaveFront:[11,0,1,""],WaveFrontAnalysis:[11,0,1,""],WavePoint:[11,0,1,""],WavePointSet:[11,0,1,""],ZernikeWaveFront:[11,0,1,""]},"optics.wavefront.Interferometer":{draw:[11,1,1,""],setTilt:[11,1,1,""],setWaveFront:[11,1,1,""]},"optics.wavefront.WaveFront":{fromFile:[11,1,1,""],getImage:[11,1,1,""],getOTF:[11,1,1,""],getPSF:[11,1,1,""],getValue:[11,1,1,""],plotImage:[11,1,1,""],plotOTF:[11,1,1,""],plotPSF:[11,1,1,""]},"optics.wavefront.WaveFrontAnalysis":{drawAberrationPlot:[11,1,1,""],fitSeidel:[11,1,1,""],fitZernike:[11,1,1,""],getWavePointSet:[11,1,1,""]},"optics.wavefront.WavePoint":{getPathLength:[11,1,1,""],getPhaseLength:[11,1,1,""],setWithRay:[11,1,1,""],setWithWaveFront:[11,1,1,""]},"optics.wavefront.WavePointSet":{add:[11,1,1,""],fitSeidel:[11,1,1,""],fitSeidelFunction:[11,1,1,""],fitZernike:[11,1,1,""],fitZernikeFunctionNine:[11,1,1,""],fitZernikeFunctionSixteen:[11,1,1,""],fitZernikeFunctionTwentyFive:[11,1,1,""],getPhaseValues:[11,1,1,""],leastSqrError:[11,1,1,""],setWithRays:[11,1,1,""],setWithWaveFront:[11,1,1,""],zeroMean:[11,1,1,""]},"optics.wavelength":{AirIndex:[12,0,1,""],Blue:[12,4,1,""],BlueColourMatch:[12,4,1,""],BlueLimit:[12,4,1,""],Cadmium_F:[12,4,1,""],CauchyIndex:[12,0,1,""],FixedIndex:[12,0,1,""],GaussianSpectrum:[12,0,1,""],Green:[12,4,1,""],GreenColourMatch:[12,4,1,""],Hydrogen_F:[12,4,1,""],InfoIndex:[12,0,1,""],MaterialIndex:[12,0,1,""],Mercury_h:[12,4,1,""],Mercury_i:[12,4,1,""],PhotopicSpectrum:[12,0,1,""],PlanckSpectrum:[12,0,1,""],Red:[12,4,1,""],RedColourMatch:[12,4,1,""],RedLimit:[12,4,1,""],RefractiveIndex:[12,0,1,""],RefractiveIndexColour:[12,2,1,""],ScotopicSpectrum:[12,0,1,""],Sellmeier:[12,0,1,""],Spectrum:[12,0,1,""],TriColourSpectrum:[12,0,1,""],WaveLength:[12,0,1,""],WavelengthColour:[12,2,1,""],getDefaultWavelength:[12,2,1,""],getDesignWavelength:[12,2,1,""],getInitialDefaultWavelength:[12,2,1,""],getInitialDesignWavelength:[12,2,1,""],setDefaultWavelength:[12,2,1,""],setDesignWavelength:[12,2,1,""],setFixedAirIndex:[12,2,1,""]},"optics.wavelength.AirIndex":{copy:[12,1,1,""]},"optics.wavelength.InfoIndex":{copy:[12,1,1,""]},"optics.wavelength.PlanckSpectrum":{setTemperature:[12,1,1,""]},"optics.wavelength.RefractiveIndex":{getNd:[12,1,1,""],getNe:[12,1,1,""],getType:[12,1,1,""],getVd:[12,1,1,""],getVe:[12,1,1,""]},"optics.wavelength.Sellmeier":{fitIndex:[12,1,1,""],getAlpha:[12,1,1,""],getLambdaZero:[12,1,1,""]},"optics.wavelength.WaveLength":{draw:[12,1,1,""],getArrayDerivatives:[12,1,1,""],getArrayValues:[12,1,1,""],getDerivative:[12,1,1,""],getValue:[12,1,1,""]},"optics.zernike":{opticalZernike:[13,2,1,""],opticalZernikeName:[13,2,1,""],radial:[13,2,1,""],zernike:[13,2,1,""]},"vector.Angle":{copy:[10,1,1,""],getDegrees:[10,1,1,""],getUnit3d:[10,1,1,""],random:[10,1,1,""],setDegrees:[10,1,1,""]},"vector.Unit3d":{copy:[10,1,1,""],getAngle:[10,1,1,""],random:[10,1,1,""],reflection:[10,1,1,""],refraction:[10,1,1,""],setPolar:[10,1,1,""],setPolarDegrees:[10,1,1,""]},"vector.Vector2d":{abs:[10,1,1,""],absCube:[10,1,1,""],absNormalised:[10,1,1,""],absSquare:[10,1,1,""],angleBetween:[10,1,1,""],copy:[10,1,1,""],distance:[10,1,1,""],distanceCube:[10,1,1,""],distanceSquare:[10,1,1,""],dot:[10,1,1,""],errorSquare:[10,1,1,""],getComplex:[10,1,1,""],inverseSquare:[10,1,1,""],isValid:[10,1,1,""],negate:[10,1,1,""],normalise:[10,1,1,""],polar:[10,1,1,""],rect:[10,1,1,""],rotate:[10,1,1,""],round:[10,1,1,""],set:[10,1,1,""],setInvalid:[10,1,1,""],setLength:[10,1,1,""]},"vector.Vector3d":{absCube:[10,1,1,""],absNormalised:[10,1,1,""],absSquare:[10,1,1,""],angleBetween:[10,1,1,""],areaBetween:[10,1,1,""],copy:[10,1,1,""],cross:[10,1,1,""],distance:[10,1,1,""],distanceCube:[10,1,1,""],distanceSquare:[10,1,1,""],dot:[10,1,1,""],errorSquare:[10,1,1,""],inverseSquare:[10,1,1,""],isValid:[10,1,1,""],negate:[10,1,1,""],normalise:[10,1,1,""],polar:[10,1,1,""],propagate:[10,1,1,""],random:[10,1,1,""],rotate:[10,1,1,""],rotateAboutX:[10,1,1,""],rotateAboutY:[10,1,1,""],rotateAboutZ:[10,1,1,""],round:[10,1,1,""],set:[10,1,1,""],setInvalid:[10,1,1,""],setLength:[10,1,1,""],setPolar:[10,1,1,""],setPolarDegrees:[10,1,1,""],unitPair:[10,1,1,""]},optics:{surface:[8,3,0,"-"],wavelength:[12,3,0,"-"]},tio:{getAngle:[9,2,1,""],getAngleDegrees:[9,2,1,""],getBool:[9,2,1,""],getComplex:[9,2,1,""],getExpandedFilename:[9,2,1,""],getFilename:[9,2,1,""],getFloat:[9,2,1,""],getInt:[9,2,1,""],getOption:[9,2,1,""],getString:[9,2,1,""],getVector2d:[9,2,1,""],getVector3d:[9,2,1,""],openFile:[9,2,1,""],setJournal:[9,2,1,""],tprint:[9,2,1,""]},vector:{Angle:[10,0,1,""],Unit3d:[10,0,1,""],Vector2d:[10,0,1,""],Vector3d:[10,0,1,""]}},objnames:{"0":["py","class","Python class"],"1":["py","method","Python method"],"2":["py","function","Python function"],"3":["py","module","Python module"],"4":["py","data","Python data"]},objtypes:{"0":"py:class","1":"py:method","2":"py:function","3":"py:module","4":"py:data"},terms:{"0000e":[0,3,6,8,11],"0mm":8,"100mm":[0,3,5],"106mm":3,"10mm":[3,5],"20mm":5,"2mm":3,"2pi":10,"30mm":5,"40mm":5,"4th":11,"55um":7,"60mm":5,"6th":11,"770l":13,"80mm":5,"8th":11,"97mm":3,"abstract":[2,7,8,11,12],"boolean":[3,9],"case":[5,9,12],"class":[1,9],"default":[0,1,2,3,5,6,7,8,9,10,11],"float":[0,2,3,4,5,6,7,8,9,10,11,12,13],"function":[1,11,13],"import":[5,7,9],"int":[0,2,3,4,7,8,9,11,12,13],"new":[2,3,5,6,7,8,10,12],"null":8,"return":[0,2,3,4,5,6,7,8,9,10,11,12,13],"switch":[7,8,9,12],"true":[3,5,6,7,9,10,11,12],"try":[3,9],"while":3,Added:10,And:9,For:[0,5],Its:9,NOT:[3,5,6,7,8,10],Not:8,The:[0,2,3,4,5,6,7,8,9,10,11,12],Then:9,There:[3,5,6,7,8,9,11,12],These:[3,4,5,8,9,10,12,13],Use:[3,6,12],Will:[2,7,10,12],__bool__:7,a_or_nd:12,abb:12,abber:11,abd:8,aberr:11,about:[3,6,8,10],abov:[5,12],abs:[8,9,10],abscub:10,abscud:10,absnormalis:10,absolut:[9,10],absolutl:9,abssquar:10,abstarct:[2,7,12],accept:[2,9],access:[5,9,10,12],account:8,achromat:3,across:[0,7,11],act:8,action:9,activ:[6,7],actual:[0,2,7,11,12],adapt:12,add:[0,3,5,6,7,11,12],addbeam:7,addcollimatedbeam:7,addcollimatedparaxialbeam:7,added:[0,3,7,9,11],addgrid:0,addit:[3,11],additi:8,addmonitor:7,addpoint:6,addsourcebeam:7,addsourceparaxialbeam:7,addtestgrid:0,advis:12,affect:8,affet:8,after:9,again:[9,12],agaist:2,aglebra:5,air:[3,12],airindex:[0,1,7],algebra:1,all:[2,3,5,7,8,9,10,12,13],allow:[5,6,7,9],alon:10,along:[0,5,7,8,11],alow:7,alpha:[6,10,12],also:[0,2,3,4,5,6,7,9,10,11,12],alter:[2,3,7,8],altern:9,alwai:[0,8,9],analy:11,analys:11,analysi:[1,5,6],ands:7,angl:[0,1,2,3,5,6,7,8,11],anglebetween:10,angular:13,ani:[5,7,8,9,11,12],annular:8,annularapertur:1,aperetur:0,apertur:[0,3,7,8,11],aperur:8,apetur:7,append:[3,9,11],appli:8,applic:[9,10],approxim:12,apretur:3,arbitrari:8,area:[6,10,11],areabetween:10,aregu:5,areprocess:[],arg:[2,3,7,9,11],argumemnt:9,argument:[2,13],aross:7,arrai:[0,7,11,12,13],asci:4,aslo:12,associ:[8,12],assoiat:7,assum:[3,5,6,7,8,9,10,11,12],ast305:12,atan2:2,attach:7,attack:7,attribut:7,augument:10,auntomat:7,auotmat:5,austin:12,autoclass:8,automat:[0,7,10,12],automormalsi:10,avaial:4,avaialbl:11,averag:[6,11],axi:[0,2,3,5,6,7,8,10,11,12],b_or_vd:12,back:[0,3,5,8,12],backfocallength:[3,5],backfocalplan:[3,5],backnodalpoint:[3,5],backpow:5,backprincipalplan:[3,5],bar:8,base:[4,5,7,8,12],basic:[1,2,3,5,11,12],beam:[0,2,7,11],been:[7,9],beep:9,befor:11,begin:3,being:[0,2,3,4,5,6,7,8,9,10,11,12],beinng:7,bell:9,belong:8,below:[3,5],bend:3,best:6,beta:10,between:[0,2,3,5,10,11],beween:[2,12],biconvex:3,bicovex:3,binari:9,bk7:[3,12],blackbodi:12,blank:[7,11],block:[7,8],blue:12,bluecolourmatch:12,bluelimit:12,bodi:12,book:11,bool:[3,5,7,9,10,11,12],born:13,both:[8,10,11],boundari:10,box:3,bright:12,bruton:12,buffer:9,bui:3,build:5,built:9,byconvex:[],cacault:12,cach:12,cacualt:6,cacul:11,cadium:12,cadmium_c:12,cadmium_f:12,cal:12,calccul:12,calcual:10,calcualt:[3,5,6,7,8,10,11,12,13],calcuat:12,calcul:[3,5,7,12,13],calculu:2,calcuuat:11,call:[3,5,7,8,9,11,12],callabl:9,calu:12,can:[3,5,7,8,9,10,11,12],cancluat:[],cannot:5,cardin:[3,5],cardinalpoint:[3,5],care:3,caucbi:12,cauchi:12,cauchyindex:1,caviti:5,cavitymatrix:5,centr:[3,6,8,11],central:3,centriod:6,centroid:6,chain:10,chang:[3,6,7,8,10,12],charact:9,check:[3,7,9,10],chevron:8,choic:9,chosen:9,circl:[6,11],circular:[0,2,3,7,8],circularapertur:[1,3,7],circularli:2,citera:3,claas:3,clase:7,clear:8,cleck:7,cloeset:0,close:9,code:[1,9,10,12],coef:[4,12],coeffic:[4,11,12],coefficeint:[4,11],coeffici:12,coefficienc:11,collim:[7,11],color:[6,12],colour:[1,6,7,11],com:12,combin:5,come:[3,7],comlour:0,command:1,comment:9,commom:5,common:[3,5,12],commonli:9,commut:5,compat:[5,8],compex:2,complex:[2,5,8,9,10,13],compoent:[2,10],compolex:[],compon:[0,2,3,5,8,10,11,13],componen:10,componet:[2,3,9,10],compont:[2,3],compound:[3,5],computation:12,comtain:9,concatin:9,consid:6,consist:[2,5,8,11],constant:[1,12],constructor:12,constuctor:11,contain:[9,11],content:1,control:[3,9,12],convent:4,convert:9,conveterd:9,convex:[],convexplano:3,coodin:5,coordin:[0,3,5,6,7,8],copi:[2,5,7,8,10,12],cordilan:3,correct:[0,8,9],coss:7,cover:12,cox:12,creat:[5,6,11],creation:[7,10],critic:10,cross:[7,10],crosseszero:7,crown:3,crownindex:3,cube:10,curatur:5,curcatur:[3,8],curfatur:8,currect:11,currectlen:3,currennt:2,current:[0,2,3,5,6,7,8,9,10,11,12],currentlen:1,curv:[5,8],curvatu:[3,5],curvatur:[3,5,8],curvatut:5,curve_fit:11,dan:12,dark:12,dat:9,data:[4,5,9,11],databas:[4,12],databaselen:1,databasematrix:5,datapoint:11,dattbas:12,deal:[9,12],dealt:7,deatil:3,decim:[9,10],deconfigur:0,deep:[7,10,12],defaault:[2,9],defaul:3,defaultnam:9,defaulttyp:9,defaultwavelength:12,defaut:[3,9],defedin:5,defein:8,defin:[3,5,7,8,10,12,13],definb:7,definit:11,defult:[2,5],degre:[9,10,12],deign:3,delta:[2,3,5,8,12],depend:[3,7,8],dependand:12,deriv:12,describ:10,descriptor:9,design:[0,1,2,3,11],designwavelength:12,det:10,detail:[5,12],detauil:3,deterim:10,determin:[2,5,8],develop:9,diagram:[1,12],diamet:3,diection:11,dielectricmatrix:5,differ:[2,12],differebnc:2,digit:12,dilectr:5,dimenensio:11,dimension:[0,3,5,6,8,10,11],dir:9,direct:[0,7,8,9,10],director:[7,8,9],dirn:7,disagram:6,disgram:6,dispers:12,displai:[0,5,9,11],distanc:[0,3,5,7,8,10],distancecub:10,distancesquar:10,divid:[],doe:[7,8,9,10],done:12,dot:10,doubl:[3,5],doublet:[1,5],doubletmatrix:5,draw:[0,3,5,6,7,8,11,12],drawaberrationplot:11,drawpsf:6,due:10,dure:7,dynam:12,each:[0,3,7,8,9,11,12],easier:[3,5],easili:9,eccenr:6,eccentr:[6,8],eddg:[],edg:[0,3,8],edgeplan:8,effic:10,effort:9,either:[3,7,8,11,12],elem:5,element:[3,5,11],elementari:9,elemnt:11,elips:6,ellip:2,ellips:[2,6],ellipt:2,els:[7,9,10,12],emissit:12,emmistt:12,empti:[0,5],encourt:3,end:3,ensur:0,entranc:3,entranceapertur:3,entrancepupil:3,env:9,environemnt:12,environment:9,environmentalvari:12,epsilon:8,equal:7,equival:[5,12],equivali:10,erasiet:[],error:[9,10,11,12],errorsquar:10,escap:1,estmat:11,etc:[7,8],etens:5,etent:0,eth:12,ether:5,evalu:[9,11],even:9,everi:[7,12],exaclti:12,exampl:[1,7,12],exceed:10,execut:9,exist:9,exit:[3,7,9,11],exitapertur:3,exitpupil:[3,11],expand:9,expans:11,expend:3,experiment:12,express:9,extarct:12,extend:[1,3,7,8,10,11,12],extens:[7,9],extent:[5,11,12],extra:[3,8,10],extract:5,eye:12,f_or_cl:5,faction:8,factor:[3,5,8,10,12,13],fail:[7,8,9,10],failur:[9,10],fals:[3,5,7,8,9,10,12],fast:2,featur:12,field:3,fig:10,figur:10,file:[1,3,5,11],filenam:[1,3,4,5],fill:[0,7],find:6,finer:9,finit:5,first:[0,3,5,6,7,11],fit:[6,11,12],fitindex:12,fitseidel:11,fitseidelfunct:11,fitzernik:11,fitzernikefunctionnin:11,fitzernikefunctionsixteen:11,fitzernikefunctiontwentyf:11,five:9,fix:[3,7,8,12],fixedairindex:12,fixedfoc:3,fixedindex:1,fixedmo:6,fixedradiu:3,flag:[3,5],flaot:10,flat:[3,5,8,11],flatsurfac:1,flintindex:3,floar:5,flocal:5,flot:[],flrat:8,flush:9,fno:3,focal:[3,5],focallength:3,follow:[5,8,9],fomula:5,forc:[10,12],form:[2,3,4,5,6,7,8,9,10,11,12,13],format:[3,5,7,9,10,12],formul:11,formula:[3,4,10,12],formult:3,fortran:12,found:5,four:2,fourier:11,franhoff:3,fred:9,fring:11,from:[0,2,3,4,5,6,7,8,9,10,11,12,13],fromat:9,fromfil:11,front:[3,5],frontfocallength:[3,5],frontfocalplan:[3,5],frontnodalpoint:[3,5],frontpow:5,frontprincipalplan:[3,5],fucntion:9,fuction:12,full:[3,7,9],fulli:2,futur:9,gamma:10,gap:0,garphic:5,gase:4,gaussian:12,gaussianspectrum:12,gener:[2,3,8,10,11,12],geometr:[3,5,6],get:[0,2,3,4,5,6,7,8,9,10,11,12],getalpha:12,getangl:[2,9,10],getangledegre:9,getarrayderiv:12,getarrayvalu:12,getbend:3,getbool:9,getcomplex:[9,10],getcurrentlen:3,getdefaultwavelength:12,getdegre:10,getderiv:12,getdesignwavelength:12,getdist:8,getedgethick:3,getellip:2,getexpandedfilenam:9,getfilenam:9,getfloat:9,getfno:3,getimag:[0,11],getinfo:[3,5,7],getinitialdefaultwavelength:12,getinitialdesignwavelength:12,getint:9,getintens:[2,7],getlambdazero:12,getlist:4,getmateri:4,getmatrix:2,getn:12,getnd:12,getnorm:8,getopt:[4,9],getotf:11,getparaxialinteract:8,getpathlength:11,getpencil:0,getphas:2,getphaselength:[7,11],getphasevalu:11,getpixelsourcepoint:0,getpoint:[3,5,8],getpsf:11,getradiu:[3,8],getraypencil:0,getsingletpair:3,getsiz:8,getsourcepoint:8,getstr:9,getsurfaceinteract:[0,8],getsystemimag:0,getter:3,getthick:3,gettyp:12,getunit3d:10,getv:12,getvalu:[11,12],getvd:12,getvector2d:9,getvector3d:9,getwavepointset:11,give:[0,3,5,8,11,12],given:[0,2,3,5,6,7,9,11,12,13],glabal:[],glass:[3,4,8],global:[0,3,5,7,8,12],globalcoordin:3,gobal:[3,8],golob:6,grai:0,graphic:1,graphm:12,green:[11,12],greencolourmatch:12,greyscal:0,grid:0,group:[3,5,7,8],group_pt:3,guassian:12,guassianspectrum:1,guess:6,gui:3,had:12,half:2,halfwav:2,halfwavepl:1,hand:[2,4,10],handel:[3,4],handl:[3,4,7,10,12],has:[3,5,7,8,9,12],have:[0,5,7,12],height:[3,5,7,8],held:[0,3,4,7,11,12],helium_d:12,hello:9,help:9,here:[2,5],hex:9,hexstr:12,high:[4,12],higher:5,histor:8,hold:[0,2,3,4,5,7,10,11,12],home:9,horizont:[0,3,8,11],hot:12,hotizont:11,href:12,html:12,http:[4,12],hydroden_c:12,hydroden_f:12,hydrodgen:12,hydrogen_f:12,idea:2,ideal:[3,11],identi:9,ignor:[0,3,11],ima:[3,5],imag:[0,3,5,8,11,12],imageplac:0,imageplan:[0,1,3,5],imagepoint:[3,5],imaginari:9,imaginig:0,imform:5,implem:2,implemend:11,implement:[0,2,3,5,6,7,8,10,11,12],implementr:5,implemnt:12,implenet:11,implet:7,implmeent:10,implment:[6,12],implwmwnt:10,imposs:[5,8],imput:9,imshow:0,in_height:5,incldu:[8,9],includ:[5,8,12],increm:5,incremend:2,increment:[2,8],incrementinputplan:5,incrementsurfac:8,indeal:2,independ:[9,12],index:[0,1,2,3,4,5,7,8,10],indic:12,induvudu:3,inf:[0,5],infinit:[0,7],info:[4,12],infoiindex:12,infoindex:[1,4],infor:5,inform:[0,8,12],inherit:7,init:[0,3,11],initialis:12,inner:8,innerradiu:8,inout:9,inplement:12,inport:5,inpout:9,input:[1,2,3,5,7,8,11],inputplan:5,inputplaneheight:5,inpyt:[],integ:12,intens:[0,1,2,6,8,11,12],intensityrai:[7,8],inter:10,interact:[0,8,9],interers:5,interfac:[3,5,12],interferomet:1,intern:[3,7,9,11,12],internet:12,interv:12,intrefac:12,intren:12,invalid:[4,7,8,10],invers:[5,10],inversesquar:10,invert:3,ion:9,iri:[3,8],irisapertur:1,irisapetur:3,isvalid:[7,10],its:[2,3,7,10],itself:7,itter:[0,3],ittera:[],ius:11,jone:1,jonematrix:2,jonesmarix:[],jonesmatix:2,jonesmatric:2,jonesmatrix:1,jonesmatrixsystem:1,jonesvector:1,jonevector:2,journal:1,just:[4,7,11,12],keep:3,kei:[0,2,3,4,7,9,11,12],kelvin:12,keyword:3,kingslak:11,kingslakewavefront:1,knife:[0,8],knifeapertur:1,knifeedgeplot:1,knifetest:0,lambda:12,lambda_0:12,larg:[11,12],laser:5,last:[3,12],law:10,lead:9,least:11,leastsqrerror:11,ledgend:[5,11],left:[2,3,5,11],leftcircularpolarisedbeam:1,legal:13,legend:[3,5,11],len:[0,1,5,8,9,11],lene:3,lenght:5,length:[3,5,7,9,10,11,12],lens:[0,3,5],lenth:5,less:[3,13],level:[5,12],lienar:3,light:[2,12],limit:12,line:[5,9,11,12],linear:[2,5],linearpolaris:1,linearpolarisedbeam:2,linearpolarsiedbeam:1,linix:9,liquid:4,list:[0,2,3,4,5,6,7,8,9,10,11,12,13],littl:5,load:4,locaion:6,local:[0,3,12],localion:0,locat:[0,3,5,6,7,8,11],locatt:3,log:11,logic:[9,11],look:12,lookup:[4,13],low:4,lower:[5,9,11],maco:9,mag:[0,3,5,10],magif:3,magificant:3,magmif:[3,5],magnif:[0,3,5],magnitud:10,mai:[3,8,9,11,12],main:[3,7,8,9,10,11,12],mainli:[3,4,12],major:6,make:[2,3,5,7,8,11,12],makestandalon:[3,8],malacara:11,mani:10,manipul:[5,11],map:0,mark:8,mask:0,mass:[3,7],match:[3,5,11,12],materail:[4,12],materi:[1,3,9,12],materialdata:[1,12],materialindex:[1,4],math:7,matplotlib:[5,6,12],matplotlinb:8,matric:[2,5],matrix:[1,2,3,7,8],matrx:[],max:[3,9],maxab:9,maximum:[5,9,11],maxium:5,maxradiu:[3,5,7,8,11],maxtrix:5,maxwavelength:12,maxwavelenth:12,maytrix:7,mean:11,measur:[5,8],mecthod:3,member:8,membership:8,menu:1,mercuri:12,mercury_:12,mercury_h:12,mercury_i:12,messag:3,methiod:[],method:[0,1,3,4,6,7,8,9,10,11,12,13],methodd:8,methoid:3,mew:5,microm:12,micron:[7,12],middl:5,miminium:9,min:9,minor:6,minotor:7,minthick:3,minumum:[],minwavelength:12,minwavelenth:12,mirror:[5,8],mirrormatrix:5,mixtur:7,mnatrix:2,modal:5,mode:0,model:12,modifi:2,modul:[1,4,9],modulu:9,moment:[1,9],mon:7,monitor:7,more:[3,5,8,9],most:[5,7,9],move:[3,8],moveabl:8,movepoint:[3,8],muiltipl:5,mult:7,multbi:7,multipl:5,multipli:[2,5,7],must:11,name:[4,9,12,13],nan:[7,10,11,13],ndarrai:[11,13],need:[7,8,11,12],neg:0,negat:[5,10],neglig:12,newlin:9,nmpy:0,nnnvvv:12,nodal:[3,5],nojourn:9,non:[3,12],nond:5,none:[0,2,3,4,5,6,7,8,9,10,11,12,13],noramlli:7,nore:11,normaal:7,normal:[0,3,5,7,8,10,12],normalis:10,normalsi:[10,11],nornmal:11,note:[0,3,4,5,6,8,9,10,11,12],noth:7,nrai:[0,7,11],number:[0,2,5,7,9,10,11,12,13],numer:[3,12],numner:[0,12],numpi:[0,10,11,12],obect:8,obj:[3,5],object:[0,3,4,5,6,7,8,12],obtain:5,oct:9,odd:[0,3],off:[7,9,11],oft:11,onc:12,one:[0,2,5,7,9,11,12],onfinit:3,onli:[0,6,7,8,9,11,12,13],oof:[],open:[3,8,9],openfil:[5,9],oper:[5,7,10],operati:[],opt:9,optain:9,optic:[0,2,3,4,5,6,7,8,9,10,11,12,13],opticalgroup:[0,1,5,7,8],opticalimag:1,opticalplan:[0,1,3,6,7],opticalsurfac:3,opticalzernik:13,opticalzernikenam:13,opticlgroup:7,optim:[6,11],optimalarea:6,option:[1,3,5,8,10,11],optitc:8,order:[2,3,5,6,7,10,11,13],org:4,orient:2,origin:[9,12],otf:11,other:[3,5,7,8,9,10,11],othwis:[],out_height:5,outer:8,outerradiu:8,output:[1,2,5,8],outputplan:5,outsid:[8,11],overal:12,overlap:5,overload:[2,5,7,10],overridd:3,overwrit:3,own:7,packag:[3,7,9,10,12],page:[1,13],pair:[0,3,5,10],pancil:7,panel:3,parabol:8,parabolicsurfac:1,paraet:5,param:[3,6,7,8,9,10,12],paramet:[0,2,3,4,5,6,7,8,9,10,11,12,13],paramt:[11,12],paramtet:3,paraxalgroup:5,paraxi:[0,1,3,5,6,11,12],paraxialdoublet:5,paraxialgroup:[1,3,7,8],paraxialmatric:5,paraxialmatrix:[1,3,7],paraxialmirror:5,paraxialplan:5,paraxialrai:1,paraxialthicklen:5,paraxialthinlen:5,parmet:12,parpag:5,pars:2,part:3,particular:12,parxial:8,parxialmatrix:7,parxialthinlen:5,pass:[0,2,7,9,11,12],path:[7,8,11],pathelength:[7,11],pathlength:[7,11],pattern:[0,11],paxaialgroup:5,paxaxialthinlen:5,paxial:3,pcov:12,peak:12,pencil:[0,6,7,11],perform:0,pertendicular:8,perzal:3,perzval:3,petzvalsum:3,phase:[2,7,11],phi:10,photop:12,photopicspectrum:1,photot:12,physic:[5,12],pick:7,piec:5,pint:9,pixel:[0,11],place:[2,3,5,7,8,10,12],plan:[5,9,11],planck:12,planckspectrum:1,plane:[0,3,5,6,7,8,10,11],planepair:[3,5],plank:12,plano:[],planoconvex:3,plate:2,plot:[2,3,5,6,7,8,11,12],plotimag:11,plototf:11,plotpoint:12,plotpsf:11,plt:[2,5,11],point:[0,1,2,3,5,7,8,10,11,12],pointinplan:7,poistion:5,polar:[2,10],polaris:2,polarplot:2,polarsi:2,polarsistaion:2,polynomi:1,popt:12,pos:[6,7,8],posit:[0,3,5,6,7,8,11],possibl:3,post:2,postit:8,potiton:5,pow:10,power:[5,7,13],pre:2,prefix:9,present:[7,8,9],princip:[3,5],principl:5,prinicp:5,prinnt:9,print:[1,5,7],printpath:7,problem:9,process:[3,7,9],processd:9,product:10,program:9,progress:7,prompt:[4,5,9],propag:[5,7,10],propagatethrough:7,propagateto:7,propagationmatrix:5,propat:7,properti:5,propgat:7,provid:9,psf:[1,11],psi:[9,10],pt_or_u:0,pt_or_x:8,pt_or_z:3,pts:6,pupil:[3,11],put:[2,5],pyplot:5,python:9,quadrat:8,quadric:8,quadricsurfac:1,quarter:2,quarterwavepl:1,qudaricsurfac:8,qudric:8,raai:11,rad:3,radi:8,radial:1,radian:[2,7,8,9,10],radiu:[0,3,5,6,7,8,11,13],radu:3,rai:[0,3,6,8,10,11],random:10,rang:[2,4,9,10,12],rather:[3,8,9],ratio:[3,8,10],ration:[8,10],raw:[],raymonitor:1,raypath:7,raypenc:0,raypencil:[0,1,6,11],raytyp:7,read:[1,3,5,11,12],readfromfil:[],readin:4,reai:10,real:[7,9,11],recalcul:12,record:[7,9],recoveri:9,rect:10,red:12,redcolourmatch:12,redlimit:12,ref:12,refactiveindex:7,refarct:[5,8,12],refarctiveindex:3,refartiveindex:3,refatc:12,refear:8,refect:10,refer:[0,3,7,8,11,12],referer:11,refix:9,reflect:[7,8,10],refopt:[0,11],reformat:9,refpt:11,refr:[7,8,12],refract:[1,3,5,7,8,10],refractic:10,refractiveindex:[3,4,7,8,12],refractiveindexcolour:12,refrat:8,refratc:[5,7,8,12],refratciveindex:[3,7],refratic:12,refrativeindex:[3,4,8,12],refrctiveindex:[],refreatc:4,refrect:12,regular:0,rel:[0,3,5,7],relat:5,releas:9,remad:3,remov:[7,8,9],removeinvalid:7,render:[6,11],repeat:5,repend:3,replac:9,repres:[0,3,5,6,8,12],represent:5,reprompt:9,request:3,requir:[5,9],rerturn:7,reset:[0,2,3,9,11,12],resolut:0,reson:12,respect:5,respes:8,respons:[9,12],rest:11,result:[3,11,12],retail:[],retain:[3,11],retand:2,retard:1,retart:2,retrurn:4,retuen:7,retuirn:[],retun:10,reurn:[],reus:9,revers:[3,8],rgb:12,right:[2,3,5,8],rightcircularpolarisedbeam:1,ring:9,root:9,rotat:[2,10],rotateabouti:10,rotateaboutx:10,rotateaboutz:10,rotatebi:2,rotatian:10,rounc:7,round:[0,10],rout:12,rubbish:10,same:[0,3,4,5,7,8,12],sampl:0,saniti:9,scalar:8,scale:[3,5,8,10,12],scipi:11,scope:9,scotop:12,scotopicspectrum:1,search:1,secifi:2,second:[3,5,6,10],see:[3,12],seidel:11,select:9,selector:9,self:[3,5,7,8,10,11,12],sellmeier:1,sellmier:12,sens:2,sensibl:[9,11],separ:5,sepcifi:3,sepectrum:12,serar:5,set:[0,2,3,5,6,7,8,9,10,11,12],setangl:2,setbend:3,setcentrethick:3,setcurrentlen:3,setcurvatur:3,setdefaultwavelength:12,setdegre:10,setdesignwavelength:12,setedgethick:3,setfixedairindex:12,setfocallength:[3,5],setfromstr:3,setinputplan:5,setinvalid:[7,10],setiri:3,setjourn:9,setknif:[0,8],setlength:10,setparamet:3,setpoint:[3,8],setpolar:10,setpolardegre:10,setradiu:3,setratio:8,setsiz:8,settemperatur:12,setthick:3,settilt:11,settin:7,setup:11,setwavefront:11,setwavelength:2,setwir:[0,8],setwithplan:[3,5],setwithrai:[6,11],setwithwavefront:11,shift:[0,2,5,8,11],shoukld:3,should:[8,9],show:[5,11],showlegend:5,side:[3,8],sign:[8,9],signific:5,simpl:[3,9,11,12],simpler:[3,8],simplesinglet:1,simplest:9,simplfi:10,simpli:[3,9,12],simul:[0,12],sinc:[3,5,9],singl:[2,3,5,7,9,11,12],singlet:1,sinmpl:6,site:12,situ:5,six:[3,12],size:[0,3,8,11],skew:7,slow:[2,11],small:[8,9],smooth:12,sopecifi:[3,5],sourc:[0,2,3,4,5,6,7,8,9,10,11,12,13],sourcebeam:7,sourceplan:7,sourcepoint:[0,1,8],space:[3,5,8,9,11],specfi:2,special:12,specid:2,specif:[0,2,5,6,7,8,9,11,12],specifi:[0,2,3,5,6,7,8,9,10,11,12,13],specifici:11,specifii:10,specral:12,spectifi:3,spectra:12,spectral:12,spectrum:[1,7],spectum:12,specturm:12,spetrum:7,sphere:10,spheric:8,sphericalsurfac:1,spot:1,spotdiagram:6,spotdigram:1,spread:1,sqr:11,squar:[0,8,10],stand:10,standalon:8,standard:[3,9,12],start:[0,3,6,7,12],startup:12,state:[2,12],stephen:12,sting:3,stop:9,str:[2,3,4,5,6,7,9,11,12,13],strike:7,string:[3,4,5,7,9],structur:9,subset:3,substract:11,success:10,sucess:7,sufac:[8,10],suitabl:7,sum:3,suppi:13,suppl:[10,12],suppli:[0,3,5,9,11],support:[1,5,9,10],sur:7,surafec:3,surfac:[0,1,3,5,7,10,11],surfaceinter:8,surfaceinteract:[0,1],surfacevector:8,surfc:7,swap:3,syntax:4,sysout:9,system:[0,2,3,5,11],take:[5,8,9,12],taken:[0,11],taret:0,target:[0,3,5,10],targetplac:0,targetplan:1,task:0,temperatur:12,tempertaur:12,tempertur:12,tenpertur:12,term:[11,12],termin:1,tessar:5,test:[0,7,8,9,10,11],text:9,tghe:5,than:[3,5,8,9],thei:[3,5,7],them:[5,10,11],theta:[2,8,9,10],thi:[0,2,3,5,6,7,8,9,10,11,12,13],thick:[0,3,5,8],thickess:3,thicklensmatrix:5,thii:8,thin:[3,5],thinest:[],thinlensmatrix:5,thinnest:3,though:5,three:[3,5,7,8,9,10,11,12],through:[0,2,5,7,10],throughout:10,throughpolaris:2,tilt:11,time:7,tio:[1,4,5],titl:[5,12],tken:[],togeth:[2,5],token:3,tprint:9,trace:[0,1,2,3,5,8,10],tragetplan:0,trail:9,transfer:8,transmiss:2,transmitt:[],trem:0,treturn:5,tri:9,triangl:10,tricolourspectrum:1,trupl:[8,9,10],tupl:[9,10],turn:7,two:[0,2,3,5,6,7,8,10,11,12],twyman:11,type:[0,1,3,4,7,8,10,11,12],typic:[3,4,5,6,7,8,9,10,12],udat:7,ued:11,unchang:5,under:[0,9],underli:[3,5,8],uniqu:[8,9],unit3d:[0,1,5,7,8,11],unit:[2,5,8,10,11],uniti:[5,10],unitpair:10,univers:12,unix:9,untest:9,upadt:7,updat:[3,7,8,12],updatemonitor:7,upper:9,urvatur:[],use:[0,3,7,8,9,10,12],used:[0,3,5,6,7,8,9,10,11,12],useful:[3,7,8,9,11,12],user:[3,5,7,9,11,12],usernam:9,uses:[3,9],using:[0,2,3,4,5,7,10,11,12],usual:[4,8,12],util:9,valid:[4,6,7,8,10,12],valu:[1,2,3,6,7,8,9,10,11,13],valueerror:8,valul:10,vari:3,variabl:[0,3,7,8,9,12],varial:12,variiabl:12,variou:[0,5,7,8,10],vaue:[],vaul:12,vecftor3d:10,vector2d:[0,1,6,7,8,9,11,13],vector3:8,vector3d:[0,1,3,5,6,7,8,9],vector:[1,2,5,7,8],vectors2d:10,vectors3d:10,vectror:10,version:8,vertic:[3,8,11],vetor3d:10,vextor3d:8,via:[0,4,5,7,8,9,11,12],view:11,virtual:11,visual:12,wafefront:11,wai:[3,8],wave:[0,2,3,6,7,11,12],wave_arrai:12,wavefront:1,wavefrontanalysi:1,wavelangth:4,wavelelength:6,wavelengh:3,wavelength:[0,1,2,3,6,7,8,11],wavelengthcolour:[7,12],wavelenth:[7,11],wavellength:2,wavelnegth:[2,3,11,12],wavelngth:3,wavelnngth:3,wavepoint:1,wavepointset:1,wavepooint:11,weavelength:2,websit:12,weight:13,what:[9,11],when:[7,9,12],where:[0,2,4,5,6,7,8,9,10,12],whhere:12,which:[3,5,7,8,9,11,12],whick:3,white:9,whith:3,whole:[3,7,12],wide:0,width:12,wih:8,wire:[0,8],within:[10,12],without:[3,8],wold:13,wolf:13,work:[5,7,9,11],world:9,wrang:[4,12],write:9,written:9,wrong:9,wrt:[2,8,10,11],wuth:8,www:12,x_or_v:10,xgap:0,xpixel:0,xsize:[0,3,8],xtilt:11,xxx:12,xxxyyi:12,yes:9,ygap:0,you:10,ypixel:0,ysize:[0,3,8],ytilt:11,zaxi:0,zerkin:13,zern:11,zernik:[1,11],zernikewavefront:1,zero:[5,9,10,11,12],zerom:[],zeromean:11,zerr:11},titles:["Analysis Classes","Welcome to Ray Optics\u2019s documentation!","Jones Methods","Lens Classes","Material Classes","Matrix Methods","Point Spread Functions and Spot Diagrams","Ray Tracing","Surface Classes","Terminal Input / Output (tio)","Vector Classes","Wavefront analysis","Wavelength Module","Zernike analysis"],titleterms:{"class":[0,2,3,4,5,6,7,8,10,11,12],"default":12,"function":[3,6,9,12],airindex:12,algebra:5,analysi:[0,11,13],angl:[9,10],annularapertur:8,basic:9,cauchyindex:12,circularapertur:8,code:5,colour:12,command:9,constant:8,currentlen:3,databaselen:3,design:12,diagram:6,document:1,doublet:3,escap:9,exampl:[5,9],extend:5,file:9,filenam:9,fixedindex:12,flatsurfac:8,graphic:5,guassianspectrum:12,halfwavepl:2,imageplan:8,index:12,indic:1,infoindex:12,input:9,intens:7,interferomet:11,irisapertur:8,jone:2,jonesmatrix:2,jonesmatrixsystem:2,jonesvector:2,journal:9,kingslakewavefront:11,knifeapertur:8,knifeedgeplot:0,leftcircularpolarisedbeam:2,len:3,linearpolaris:2,linearpolarsiedbeam:2,materi:4,materialdata:4,materialindex:12,matrix:5,menu:9,method:[2,5],modul:12,moment:6,optic:1,opticalgroup:3,opticalimag:0,opticalplan:8,option:9,output:9,parabolicsurfac:8,paraxi:7,paraxialgroup:5,paraxialmatrix:5,paraxialrai:7,photopicspectrum:12,planckspectrum:12,point:6,polynomi:13,print:9,psf:6,quadricsurfac:8,quarterwavepl:2,radial:13,rai:[1,7],raymonitor:7,raypencil:7,read:9,refract:12,retard:2,rightcircularpolarisedbeam:2,scotopicspectrum:12,sellmeier:12,simplesinglet:3,singlet:3,sourcepoint:7,spectrum:12,sphericalsurfac:8,spot:6,spotdigram:6,spread:6,support:12,surfac:8,surfaceinteract:8,tabl:1,targetplan:0,termin:9,tio:9,trace:7,tricolourspectrum:12,type:9,unit3d:10,valu:12,vector2d:10,vector3d:10,vector:[9,10],wavefront:11,wavefrontanalysi:11,wavelength:12,wavepoint:11,wavepointset:11,welcom:1,zernik:13,zernikewavefront:11}})