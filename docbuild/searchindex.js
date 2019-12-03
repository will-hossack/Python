Search.setIndex({docnames:["analysis","index","lens","matrix","psf","ray","surface","tio","vector","wavelength","zernike"],envversion:{"sphinx.domains.c":1,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":1,"sphinx.domains.javascript":1,"sphinx.domains.math":2,"sphinx.domains.python":1,"sphinx.domains.rst":1,"sphinx.domains.std":1,"sphinx.ext.intersphinx":1,"sphinx.ext.todo":2,"sphinx.ext.viewcode":1,sphinx:56},filenames:["analysis.rst","index.rst","lens.rst","matrix.rst","psf.rst","ray.rst","surface.rst","tio.rst","vector.rst","wavelength.rst","zernike.rst"],objects:{"optics.analysis":{AberrationPlot:[0,0,1,""],KnifeEdgeTest:[0,0,1,""],OpticalImage:[0,0,1,""],TargetPlane:[0,0,1,""],WavePoint:[0,0,1,""],WavePointSet:[0,0,1,""]},"optics.analysis.AberrationPlot":{draw:[0,1,1,""]},"optics.analysis.KnifeEdgeTest":{getImage:[0,1,1,""],setKnife:[0,1,1,""]},"optics.analysis.OpticalImage":{addTestGrid:[0,1,1,""],draw:[0,1,1,""],getImage:[0,1,1,""],getPixelSourcePoint:[0,1,1,""],getRayPencil:[0,1,1,""],getSurfaceInteraction:[0,1,1,""],getSystemImage:[0,1,1,""]},"optics.analysis.TargetPlane":{add:[0,1,1,""],addGrid:[0,1,1,""],draw:[0,1,1,""],getPencils:[0,1,1,""],rayPencil:[0,1,1,""]},"optics.analysis.WavePoint":{getPhaseLength:[0,1,1,""],setWithRay:[0,1,1,""],setWithZernike:[0,1,1,""]},"optics.analysis.WavePointSet":{add:[0,1,1,""],fitZernike:[0,1,1,""],getPhaseValues:[0,1,1,""],leastSqrError:[0,1,1,""],setWithPencil:[0,1,1,""],setWithZernike:[0,1,1,""],zeroMean:[0,1,1,""]},"optics.lens":{DataBaseLens:[2,0,1,""],Doublet:[2,0,1,""],Lens:[2,0,1,""],OpticalGroup:[2,0,1,""],SimpleSinglet:[2,0,1,""],Singlet:[2,0,1,""]},"optics.lens.Doublet":{getFNo:[2,1,1,""],getRadius:[2,1,1,""],invert:[2,1,1,""],setCurvatures:[2,1,1,""],setFocalLength:[2,1,1,""],setRadius:[2,1,1,""]},"optics.lens.Lens":{backFocalLength:[2,1,1,""],backFocalPlane:[2,1,1,""],backNodalPoint:[2,1,1,""],backPrincipalPlane:[2,1,1,""],cardinalPoints:[2,1,1,""],draw:[2,1,1,""],entrancePupil:[2,1,1,""],exitPupil:[2,1,1,""],frontFocalLength:[2,1,1,""],frontFocalPlane:[2,1,1,""],frontNodalPoint:[2,1,1,""],frontPrincipalPlane:[2,1,1,""],petzvalSum:[2,1,1,""],setFocalLength:[2,1,1,""],setIris:[2,1,1,""]},"optics.lens.OpticalGroup":{add:[2,1,1,""],draw:[2,1,1,""],entranceAperture:[2,1,1,""],exitAperture:[2,1,1,""],getInfo:[2,1,1,""],imagePoint:[2,1,1,""],paraxialGroup:[2,1,1,""],paraxialMatrix:[2,1,1,""],planePair:[2,1,1,""],scale:[2,1,1,""],setPoint:[2,1,1,""]},"optics.lens.Singlet":{getBend:[2,1,1,""],getEdgeThickness:[2,1,1,""],getFNo:[2,1,1,""],getRadius:[2,1,1,""],getThickness:[2,1,1,""],setBend:[2,1,1,""],setCentreThickness:[2,1,1,""],setCurvatures:[2,1,1,""],setEdgeThickness:[2,1,1,""],setFocalLength:[2,1,1,""],setFromString:[2,1,1,""],setParameters:[2,1,1,""],setRadius:[2,1,1,""],setThickness:[2,1,1,""]},"optics.matrix":{CavityMatrix:[3,0,1,""],DataBaseMatrix:[3,0,1,""],DielectricMatrix:[3,0,1,""],DoubletMatrix:[3,0,1,""],MirrorMatrix:[3,0,1,""],ParaxialDoublet:[3,0,1,""],ParaxialGroup:[3,0,1,""],ParaxialMatrix:[3,0,1,""],ParaxialMirror:[3,0,1,""],ParaxialPlane:[3,0,1,""],ParaxialThickLens:[3,0,1,""],ParaxialThinLens:[3,0,1,""],PropagationMatrix:[3,0,1,""],ThickLensMatrix:[3,0,1,""],ThinLensMatrix:[3,0,1,""]},"optics.matrix.ParaxialGroup":{backFocalPlane:[3,1,1,""],backNodalPoint:[3,1,1,""],backPrincipalPlane:[3,1,1,""],cardinalPoints:[3,1,1,""],copy:[3,1,1,""],draw:[3,1,1,""],frontFocalPlane:[3,1,1,""],frontNodalPoint:[3,1,1,""],frontPrincipalPlane:[3,1,1,""],getInfo:[3,1,1,""],getPoint:[3,1,1,""],imagePlane:[3,1,1,""],imagePoint:[3,1,1,""],inputPlane:[3,1,1,""],maxRadius:[3,1,1,""],outputPlane:[3,1,1,""],planePair:[3,1,1,""],scale:[3,1,1,""]},"optics.matrix.ParaxialMatrix":{backFocalLength:[3,1,1,""],backFocalPlane:[3,1,1,""],backPower:[3,1,1,""],backPrincipalPlane:[3,1,1,""],copy:[3,1,1,""],determinant:[3,1,1,""],frontFocalLength:[3,1,1,""],frontFocalPlane:[3,1,1,""],frontPower:[3,1,1,""],frontPrincipalPlane:[3,1,1,""],inverse:[3,1,1,""],scale:[3,1,1,""],setFocalLength:[3,1,1,""],trace:[3,1,1,""]},"optics.matrix.ParaxialPlane":{draw:[3,1,1,""],getInfo:[3,1,1,""]},"optics.psf":{FixedMoments:[4,0,1,""],Psf:[4,0,1,""],SpotDiagram:[4,0,1,""]},"optics.psf.FixedMoments":{addPoint:[4,1,1,""],area:[4,1,1,""],centroid:[4,1,1,""],eccentricity:[4,1,1,""],ellipse:[4,1,1,""],radius:[4,1,1,""]},"optics.psf.Psf":{area:[4,1,1,""],draw:[4,1,1,""],eccentricity:[4,1,1,""],ellipse:[4,1,1,""],optimalArea:[4,1,1,""],setWithRays:[4,1,1,""]},"optics.psf.SpotDiagram":{draw:[4,1,1,""]},"optics.ray":{IntensityRay:[5,0,1,""],ParaxialRay:[5,0,1,""],PrintPath:[5,0,1,""],Ray:[5,0,1,""],RayMonitor:[5,0,1,""],RayPath:[5,0,1,""],RayPencil:[5,0,1,""],SourcePoint:[5,0,1,""]},"optics.ray.IntensityRay":{copy:[5,1,1,""],getPhaselength:[5,1,1,""],isValid:[5,1,1,""],pointInPlane:[5,1,1,""],propagate:[5,1,1,""],propagateThrough:[5,1,1,""],setInvalid:[5,1,1,""]},"optics.ray.ParaxialRay":{copy:[5,1,1,""],crosses:[5,1,1,""],crossesZero:[5,1,1,""],isValid:[5,1,1,""],mult:[5,1,1,""],multBy:[5,1,1,""],propagate:[5,1,1,""],propagateThrough:[5,1,1,""],propagateTo:[5,1,1,""],setInvalid:[5,1,1,""]},"optics.ray.PrintPath":{update:[5,1,1,""]},"optics.ray.Ray":{addMonitor:[5,1,1,""],draw:[5,1,1,""],isValid:[5,1,1,""],setInvalid:[5,1,1,""],updateMonitor:[5,1,1,""]},"optics.ray.RayPath":{copy:[5,1,1,""],draw:[5,1,1,""],getInfo:[5,1,1,""],update:[5,1,1,""]},"optics.ray.RayPencil":{addCollimatedBeam:[5,1,1,""],addCollimatedParaxialBeam:[5,1,1,""],addMonitor:[5,1,1,""],addSourceBeam:[5,1,1,""],addSourceParaxialBeam:[5,1,1,""],draw:[5,1,1,""],propagate:[5,1,1,""],propagateThrough:[5,1,1,""],removeInvalid:[5,1,1,""]},"optics.ray.SourcePoint":{clone:[5,1,1,""],copy:[5,1,1,""],getIntensity:[5,1,1,""]},"optics.surface":{AnnularAperture:[6,0,1,""],CircularAperture:[6,0,1,""],Clear:[6,3,1,""],FlatSurface:[6,0,1,""],ImagePlane:[6,0,1,""],IrisAperture:[6,0,1,""],KnifeEdgeAperture:[6,0,1,""],OpticalPlane:[6,0,1,""],ParabolicSurface:[6,0,1,""],QuadricSurface:[6,0,1,""],Reflecting:[6,3,1,""],Refracting:[6,3,1,""],SphericalSurface:[6,0,1,""],Surface:[6,0,1,""]},"optics.surface.AnnularAperture":{draw:[6,1,1,""],getNormal:[6,1,1,""],getSurfaceInteraction:[6,1,1,""],scale:[6,1,1,""]},"optics.surface.CircularAperture":{draw:[6,1,1,""],getNormal:[6,1,1,""],getSurfaceInteraction:[6,1,1,""],scale:[6,1,1,""]},"optics.surface.FlatSurface":{getDistance:[6,1,1,""],getNormal:[6,1,1,""],getSurfaceInteraction:[6,1,1,""]},"optics.surface.ImagePlane":{draw:[6,1,1,""],scale:[6,1,1,""],setSize:[6,1,1,""]},"optics.surface.IrisAperture":{draw:[6,1,1,""],getNormal:[6,1,1,""],getSurfaceInteraction:[6,1,1,""],setRatio:[6,1,1,""]},"optics.surface.KnifeEdgeAperture":{getSurfaceInteraction:[6,1,1,""],setKnife:[6,1,1,""]},"optics.surface.OpticalPlane":{draw:[6,1,1,""],getSourcePoint:[6,1,1,""],getSurfaceInteraction:[6,1,1,""],incrementSurface:[6,1,1,""],surfaceVector:[6,1,1,""]},"optics.surface.QuadricSurface":{draw:[6,1,1,""],eccentricity:[6,1,1,""],edgePlane:[6,1,1,""],getDistance:[6,1,1,""],getNormal:[6,1,1,""],getParaxialInteraction:[6,1,1,""],getSurfaceInteraction:[6,1,1,""],scale:[6,1,1,""]},"optics.surface.Surface":{draw:[6,1,1,""],getNormal:[6,1,1,""],getPoint:[6,1,1,""],getSurfaceInteraction:[6,1,1,""],makeStandAlone:[6,1,1,""],scale:[6,1,1,""],setPoint:[6,1,1,""]},"optics.wavelength":{AirIndex:[9,0,1,""],Blue:[9,3,1,""],BlueColourMatch:[9,3,1,""],BlueLimit:[9,3,1,""],Cadmium_F:[9,3,1,""],CauchyIndex:[9,0,1,""],FixedIndex:[9,0,1,""],GaussianSpectrum:[9,0,1,""],Green:[9,3,1,""],GreenColourMatch:[9,3,1,""],Hydrogen_F:[9,3,1,""],InfoIndex:[9,0,1,""],MaterialIndex:[9,0,1,""],Mercury_h:[9,3,1,""],Mercury_i:[9,3,1,""],PhotopicSpectrum:[9,0,1,""],PlanckSpectrum:[9,0,1,""],Red:[9,3,1,""],RedColourMatch:[9,3,1,""],RedLimit:[9,3,1,""],RefractiveIndex:[9,0,1,""],RefractiveIndexColour:[9,4,1,""],ScotopicSpectrum:[9,0,1,""],Spectrum:[9,0,1,""],TriColourSpectrum:[9,0,1,""],WaveLength:[9,0,1,""],WavelengthColour:[9,4,1,""],getDefaultWavelength:[9,4,1,""],setDefaultWavelength:[9,4,1,""],setFixedAirIndex:[9,4,1,""]},"optics.wavelength.AirIndex":{copy:[9,1,1,""]},"optics.wavelength.InfoIndex":{copy:[9,1,1,""]},"optics.wavelength.PlanckSpectrum":{setTemperature:[9,1,1,""]},"optics.wavelength.RefractiveIndex":{getNd:[9,1,1,""],getNe:[9,1,1,""],getType:[9,1,1,""],getVd:[9,1,1,""],getVe:[9,1,1,""]},"optics.wavelength.WaveLength":{draw:[9,1,1,""],getArrayDerivatives:[9,1,1,""],getArrayValues:[9,1,1,""],getDerivative:[9,1,1,""],getValue:[9,1,1,""]},"optics.zernike":{ZernikeExpansion:[10,0,1,""],opticalZernike:[10,4,1,""],radial:[10,4,1,""],zernike:[10,4,1,""]},"optics.zernike.ZernikeExpansion":{draw:[10,1,1,""],getImage:[10,1,1,""],getValue:[10,1,1,""]},"vector.Angle":{copy:[8,1,1,""],getDegrees:[8,1,1,""],getUnit3d:[8,1,1,""],random:[8,1,1,""],setDegrees:[8,1,1,""]},"vector.Unit3d":{copy:[8,1,1,""],getAngle:[8,1,1,""],random:[8,1,1,""],reflection:[8,1,1,""],refraction:[8,1,1,""],setPolar:[8,1,1,""]},"vector.Vector2d":{abs:[8,1,1,""],absCube:[8,1,1,""],absNormalised:[8,1,1,""],absSquare:[8,1,1,""],angleBetween:[8,1,1,""],copy:[8,1,1,""],distance:[8,1,1,""],distanceCube:[8,1,1,""],distanceSquare:[8,1,1,""],dot:[8,1,1,""],errorSquare:[8,1,1,""],getComplex:[8,1,1,""],inverseSquare:[8,1,1,""],isValid:[8,1,1,""],negate:[8,1,1,""],normalise:[8,1,1,""],polar:[8,1,1,""],rect:[8,1,1,""],rotate:[8,1,1,""],round:[8,1,1,""],set:[8,1,1,""],setInvalid:[8,1,1,""],setLength:[8,1,1,""]},"vector.Vector3d":{absCube:[8,1,1,""],absNormalised:[8,1,1,""],absSquare:[8,1,1,""],angleBetween:[8,1,1,""],areaBetween:[8,1,1,""],copy:[8,1,1,""],cross:[8,1,1,""],distance:[8,1,1,""],distanceCube:[8,1,1,""],distanceSquare:[8,1,1,""],dot:[8,1,1,""],errorSquare:[8,1,1,""],inverseSquare:[8,1,1,""],isValid:[8,1,1,""],negate:[8,1,1,""],normalise:[8,1,1,""],polar:[8,1,1,""],propagate:[8,1,1,""],random:[8,1,1,""],rotate:[8,1,1,""],rotateAboutX:[8,1,1,""],rotateAboutY:[8,1,1,""],rotateAboutZ:[8,1,1,""],round:[8,1,1,""],set:[8,1,1,""],setInvalid:[8,1,1,""],setLength:[8,1,1,""],setPolar:[8,1,1,""],setPolarDegrees:[8,1,1,""],unitPair:[8,1,1,""]},optics:{surface:[6,2,0,"-"],wavelength:[9,2,0,"-"]},tio:{getAngle:[7,4,1,""],getBool:[7,4,1,""],getComplex:[7,4,1,""],getExpandedFilename:[7,4,1,""],getFilename:[7,4,1,""],getFloat:[7,4,1,""],getInt:[7,4,1,""],getOption:[7,4,1,""],getString:[7,4,1,""],getVector2d:[7,4,1,""],getVector3d:[7,4,1,""],openFile:[7,4,1,""],setJournal:[7,4,1,""],tprint:[7,4,1,""]},vector:{Angle:[8,0,1,""],Unit3d:[8,0,1,""],Vector2d:[8,0,1,""],Vector3d:[8,0,1,""]}},objnames:{"0":["py","class","Python class"],"1":["py","method","Python method"],"2":["py","module","Python module"],"3":["py","data","Python data"],"4":["py","function","Python function"]},objtypes:{"0":"py:class","1":"py:method","2":"py:module","3":"py:data","4":"py:function"},terms:{"0000e":[0,2,4,6],"100mm":3,"10mm":3,"20mm":3,"2pi":8,"30mm":3,"36mm":0,"55um":5,"abstract":[5,6,9],"boolean":7,"case":7,"class":[1,7],"default":[0,1,2,3,4,5,6,7,8,10],"float":[0,2,3,4,5,6,7,8,9,10],"function":[1,10],"import":[3,5,7],"int":[0,2,5,6,7,9,10],"new":[3,4,5,6,8,9],"null":6,"return":[0,2,3,4,5,6,7,8,9,10],"switch":[5,7,9],"true":[2,4,5,7,8,9],"try":7,Added:8,And:7,For:[0,3],Its:7,NOT:[2,3,4,5,6,8,9],Not:6,The:[0,2,3,4,5,6,7,8,9,10],Then:7,There:[3,4,5,7,9],These:[7,8],Use:[2,4],Will:[5,8,9],__bool__:5,a_or_nd:9,abb:9,abd:6,aberr:1,aberrationplot:0,about:[4,8],abov:[3,9],abs:[6,7,8],abscub:8,abscud:8,absnormalis:8,absolut:[7,8],absolutl:7,abssquar:8,abstarct:[5,9],accept:7,access:[3,7,8],account:6,achromat:2,across:[0,5],action:7,activ:[4,5],actual:[0,5,9],adapt:9,add:[0,2,3,4,5,9],addcollimatedbeam:5,addcollimatedparaxialbeam:5,added:[0,2,5,7],addgrid:0,addit:[0,2],additi:6,addmonitor:5,addpoint:4,addsourcebeam:5,addsourceparaxialbeam:5,addtestgrid:0,affect:6,affet:6,after:7,again:[7,9],aglebra:3,air:[2,9],airindex:[0,1,5],algebra:1,all:[2,5,6,7,8,9],allow:[3,4,5,7],along:[3,5,6],alow:5,alpha:[4,8],also:[0,2,3,4,5,7],alter:[2,5,6],altern:7,alwai:[0,6,7],analys:0,analysi:[1,3,4],ands:5,angl:[0,1,2,3,4,5,6],anglebetween:8,angular:10,ani:[3,5,6,7,9],anlong:6,annular:6,annularapertur:1,aperetur:0,apertur:[0,2,5,6],aperur:6,append:[0,2,7],appli:[0,6],applic:7,approxim:9,apretur:2,arbitrari:6,area:[4,8],areabetween:8,aregu:3,areprocess:2,arg:[0,2,5,7,10],argumemnt:7,argument:10,aross:5,arrai:[0,5,9,10],aslo:9,associ:9,assoiat:5,assum:[0,2,3,4,5,6,7,8,9],ast305:9,attach:5,attack:5,attch:5,attribut:5,augument:8,auntomat:5,austin:9,automat:[5,8,9],automatical:0,automormalsi:8,averag:[0,4],axi:[0,2,3,4,5,6,8,9],b_or_vd:9,back:[0,2,3,6],backfocallength:[2,3],backfocalplan:[2,3],backnodalpoint:[2,3],backpow:3,backprincipalplan:[2,3],bar:6,base:[3,5,6,9],basic:[1,2,3,9],beam:5,been:[5,7],beep:7,begin:2,being:[0,2,3,4,5,7,8,9],beinng:5,bell:7,belong:6,below:2,bend:2,best:4,beta:8,between:[0,2,8],beween:9,biconvex:2,binari:7,bk7:[2,9],blackbodi:9,blank:10,block:[5,6],blue:9,bluecolourmatch:9,bluelimit:9,bodi:9,bool:[3,5,7,8,9],born:10,both:[6,8,10],boundari:8,bright:9,bruton:9,buffer:7,bui:2,build:3,built:7,byconvex:2,cacault:9,cach:9,cacualt:4,cacul:0,cadium:9,cadmium_c:9,cadmium_f:9,cal:9,calccul:9,calcual:8,calcualt:[0,2,3,4,5,8,9,10],calcuat:9,calcul:[2,3,5,9,10],call:[0,2,3,5,6,7,9],callabl:7,calu:9,can:[2,3,5,6,7,8,9],cannot:3,cardin:[2,3],cardinalpoint:[2,3],care:2,caucbi:9,cauchi:9,cauchyindex:1,caviti:3,cavitymatrix:3,centr:[2,4,6],centriod:4,centroid:4,chain:8,chang:[2,4,5,8,9],charact:7,check:[5,7,8],chevron:6,choic:7,chosen:7,circl:4,circular:[0,2,5,6],circularapertur:[1,2,5],claas:2,clase:5,clear:6,cleck:5,cloeset:0,clone:5,close:7,code:[1,7,9],coef:9,coeffic:9,coeffici:9,coefficienc:10,collim:5,color:[4,9],colour:[0,1,4,5],com:9,combin:3,comlour:0,command:1,comment:7,commom:3,common:[2,3,9],commonli:7,commut:3,compat:3,complex:[3,6,7,8,10],compoent:[8,10],compon:[0,2,3,6,8],componen:8,componet:[7,8],compound:2,computation:9,comtain:7,concatin:7,consid:4,consist:[3,6],constant:[1,9],constructor:9,contain:[0,2,7],content:1,control:[7,9],convert:7,conveterd:7,convex:2,convexplano:2,coodin:3,coordin:[0,2,3,4,5,6],copi:[3,5,8,9],cordilan:2,correct:[0,6,7],coss:5,cover:9,cox:9,creat:[0,3,4,9],creation:[5,8],critic:8,cross:[5,8],crosseszero:5,crownindex:2,cube:8,curatur:3,curcatur:[2,6],curfatur:6,currect:0,current:[0,2,3,4,5,6,7,8,9],curv:[3,6],curvatu:[2,3],curvatur:[2,3,6],curvatut:3,dan:9,dark:9,dat:7,data:[0,3,7,10],databas:9,databaselen:1,databasematrix:3,datapoint:0,dattbas:9,deal:[7,9],dealt:5,deatil:2,decim:[7,8],deconfigur:0,deep:[5,8,9],defaault:7,defaultnam:7,defaulttyp:7,defaultwavelength:9,defaut:2,defauttyp:7,defedin:3,defein:6,defin:[2,3,5,6,8,9,10],definb:5,defult:3,degre:[8,9],delta:[6,9],depend:[2,5,6],dependand:9,deriv:9,descriptor:7,design:0,det:8,detail:[3,9],detauil:2,deterim:8,determin:[3,6],develop:7,diagram:[1,9],diamet:2,dielectricmatrix:3,differ:9,digit:9,dilectr:3,dimension:[0,2,3,4,6,8],dir:7,direct:[0,5,6,7],director:[0,2,5,6,7],dirn:5,disagram:4,disgram:4,dispers:9,displai:[0,3,7],distanc:[0,3,5,6,8],distancecub:8,distancesquar:8,doe:[5,6,8],dot:8,doubl:[2,3],doublet:[1,3],doubletmatrix:3,draw:[0,2,3,4,5,6,9,10],drawpsf:4,due:8,dure:5,dynam:9,each:[0,2,5,7,9],eack:0,easier:2,easili:7,eccenr:4,eccentr:[4,6],eddg:2,edeg:0,edg:[0,2,6],edgeplan:6,effic:8,effort:7,either:[2,6,9],elem:3,element:[2,3],elementari:7,elemnt:0,elips:4,ellips:4,els:[5,7,8,9],emissit:9,emmistt:9,empti:[0,3],encourt:2,end:2,ensur:0,entranc:[0,2],entranceapertur:2,entrancepupil:2,env:7,environemnt:9,environment:7,environmentalvari:9,epsilon:6,equal:5,equival:3,equivali:8,erasiet:3,error:[0,7,8,9],errorsquar:8,escap:1,estmat:0,etens:3,etent:0,eth:9,ether:3,etx:5,evalu:7,even:7,everi:[5,9],exampl:[1,5,9],exceed:8,execut:7,exist:7,exit:[2,5,7],exitapertur:2,exitpupil:2,expand:7,expans:[0,10],expend:2,express:7,extarct:9,extend:[0,1,2,5,6,8,9],extens:[5,7],extent:[3,9,10],extra:[2,6,8],extract:3,eye:9,f_or_cl:3,faction:6,factor:[2,3,6,8,9,10],fail:[5,6,7,8],failur:[7,8],fals:[3,5,7,8,9],featur:9,field:2,fig:8,figur:8,file:[1,2,3],filenam:[1,3],fill:[0,5],film:6,find:4,finer:7,finit:3,first:[0,2,3,4,5],fit:[0,4],fitzernik:0,five:7,fix:[6,9],fixedairindex:9,fixedindex:1,fixedmo:4,flag:3,flaot:8,flat:[0,3,6],flatsurfac:1,flintindex:2,floar:[3,7],flocal:3,flot:2,flrat:6,flush:7,fno:2,focal:[2,3],focallength:2,follow:[3,6,7],fomula:3,forc:[8,9],form:[0,2,3,4,5,6,7,8,9,10],format:[2,3,5,7,8,9],formula:[2,8,9],fortran:9,found:3,fred:7,from:[0,2,3,4,5,6,7,8,9],fromat:7,front:[2,3],frontfocallength:[2,3],frontfocalplan:[2,3],frontnodalpoint:[2,3],frontpow:3,frontprincipalplan:[2,3],fucntion:7,full:[2,5,7],futur:7,gamma:8,gap:0,gaussian:9,gaussianspectrum:9,gener:[2,6,8,9],geometr:[2,3,4],get:[0,2,3,4,5,6,7,8,9,10],getangl:[7,8],getarrayderiv:9,getarrayvalu:9,getbend:2,getbool:7,getcomplex:[7,8],getdefaultwavelength:9,getdegre:8,getderiv:9,getdist:6,getedgethick:2,getexpandedfilenam:7,getfilenam:7,getfloat:7,getfno:2,getimag:[0,10],getinfo:[2,3,5],getint:7,getintens:5,getn:9,getnd:9,getnorm:6,getopt:7,getparaxialinteract:6,getpencil:0,getphaselength:[0,5],getphasevalu:0,getpixelsourcepoint:0,getpoint:[3,6],getradiu:2,getraypencil:0,getsourcepoint:6,getstr:7,getsurfaceinteract:[0,6],getsystemimag:0,getthick:2,gettyp:9,getunit3d:8,getv:9,getvalu:[9,10],getvd:9,getvector2d:7,getvector3d:7,give:[0,2,3,6,9],given:[0,2,3,4,5,7,9,10],glabal:3,glass:2,global:[0,2,3,5,6,9],globalcoordin:2,globar:9,gobal:[2,6],golob:4,grai:0,graphm:9,green:9,greencolourmatch:9,greyscal:0,grid:0,group:[2,3,5,6],group_pt:2,guassian:9,guassianspectrum:1,guess:4,handel:2,handl:[2,5,8,9],has:[2,5,6,7,9],have:[0,3,5,9],height:[2,3,5,6],held:[0,5],helium_d:9,hello:7,help:7,hex:7,hexstr:9,high:9,histor:6,hold:[0,2,3,5,8,9,10],home:7,horizont:[0,6],hot:9,href:9,html:9,http:9,hydroden_c:9,hydroden_f:9,hydrodgen:9,hydrogen_f:9,ideal:2,identi:7,ignor:[0,2],ima:3,imag:[0,2,3,6,9,10],imageplac:0,imageplan:[0,1,3],imagepoint:[2,3],imaginari:7,imaginig:0,imform:3,implement:[0,2,3,4,5,6,8,9],implementr:3,implemnt:[8,9],implenet:0,implet:5,implmeent:8,implment:[4,9],imposs:6,imput:7,imshow:0,in_height:3,incldu:7,includ:[3,9],increment:6,incrementsurfac:6,independ:[7,9],index:[0,1,2,3,5,6,8],indic:9,induvudu:2,inf:[0,3],infinit:[0,5],info:9,infoindex:1,infor:3,inform:[0,6,9],init:0,initialis:9,inner:6,innerradiu:6,inout:7,inport:3,inpout:7,input:[1,2,3,5,6],inputplan:3,inputplaneheight:3,integ:9,intens:[0,1,4,6,9],intensityrai:[5,6],interact:[0,6,7],interers:3,interfac:[2,3,9],interferomet:10,intern:[2,5,7,9],internet:9,interv:9,intesn:0,intrefac:9,invalid:[5,8],invers:[3,8],inversesquar:8,invert:2,ion:7,iri:[2,6],irisapertur:1,irisapetur:2,isvalid:[5,8],its:[2,5,8],itself:5,itter:0,ittera:2,ius:0,journal:1,just:[5,9],keep:2,kei:[0,2,5,7,9],kelvin:9,keyword:2,knife:[0,6],knifeedgeapertur:1,knifeedgeplot:1,knifeedgetest:0,lambda:9,larg:9,laser:3,last:[2,9],law:8,lead:7,least:0,leastsqrerror:0,left:[2,3],legal:10,legend:3,len:[0,1,3,7],lene:2,length:[0,2,3,5,7,8],lens:[0,2,3],lenth:3,less:2,level:9,lienar:2,light:9,limit:9,line:[0,7,9],linear:3,linix:7,list:[0,2,3,4,5,6,7,8,9,10],locaion:4,local:[0,9],localion:0,locat:[0,2,3,4,5,6,10],logic:7,look:9,lower:7,maco:7,mag:[0,2,3,8],magif:2,magnif:[0,3],magnitud:8,mai:[0,2,6,7,9,10],main:[2,5,7,9],mainli:9,major:4,make:[2,3,5,6,9],makestandalon:[2,6],manipul:[3,8],map:0,mark:6,mask:0,mass:[2,5],match:9,materail:9,materi:[2,7,9],materialdata:9,materialindex:1,math:5,matplotlib:[4,9],matplotlinb:6,matric:3,matrix:[1,2,5],matrx:3,max:[2,7],maxab:7,maximum:[3,7],maxium:3,maxradiu:[0,3,5,6],maxtrix:3,maxwavelength:9,maxwavelenth:9,maytrix:5,mean:0,measur:6,mecthod:2,membership:6,menu:1,mercuri:9,mercury_:9,mercury_h:9,mercury_i:9,messag:2,method:[0,1,2,4,5,6,7,8,9,10],methodd:6,mew:3,microm:9,micron:[5,9],middl:3,miminium:7,min:7,minor:4,minotor:5,minumum:2,minwavelength:9,minwavelenth:9,mirror:3,mirrormatrix:3,mixtur:5,modal:3,model:9,modul:[1,7],modulu:7,moment:[1,7],mon:5,monitor:5,more:7,most:[3,5,7],move:2,moveabl:6,muiltipl:3,mult:5,multbi:5,multipl:3,multipli:[3,5],name:[7,9],nan:[5,8,10],ndarrai:10,need:[0,5,6,9],neg:0,negat:[3,8],neglig:9,newlin:7,nmpy:0,nnnvvv:9,nodal:[2,3],nojourn:7,non:[2,9],nond:3,none:[0,2,3,4,5,6,7,8,9,10],noramlli:5,normaal:5,normal:[0,2,3,5,6,8,9],normalis:8,normalsi:8,note:[0,2,3,4,6,7,8,9],noth:5,nrai:[0,5],number:[0,3,5,7,8,9],numer:[2,9],numner:[0,9],numpi:[0,9],obj:3,object:[0,2,3,4,5,9],obtain:3,oct:7,odd:[0,2],off:[0,5,7],onc:9,one:[0,3,5,7,9],onfinit:2,onli:[0,3,4,6,7,9,10],open:[6,7],openfil:[3,7],oper:[3,5],operati:3,opt:7,optain:7,optic:[0,2,3,4,5,6,7,8,9,10],opticalgroup:[0,1,3,5,6],opticalimag:1,opticalplan:[0,1,2,4,5],opticalsurfac:2,opticalzernik:10,opticlgroup:5,optim:[0,4],optimalarea:4,option:[0,1,8],optitc:6,order:[2,3,4,5,8,10],origin:7,oth:5,other:[2,3,5,7],othwis:10,out_height:3,outer:6,outerradiu:6,output:[1,3,6,10],outputplan:3,outsid:6,overal:9,overload:3,own:5,packag:[5,7,9],page:[1,10],pair:[0,2,3,8],pancil:5,parabol:6,parabolicsurfac:1,paraet:3,param:[0,2,4,5,6,7,8,9],paramet:[0,2,3,4,5,6,7,8,9,10],paramt:9,paramtet:2,paraxi:[0,1,2,3,4],paraxialdoublet:3,paraxialgroup:[1,2,5],paraxialmatric:3,paraxialmatrix:[1,2,5],paraxialmirror:3,paraxialplan:3,paraxialrai:1,paraxialthicklen:3,paraxialthinlen:3,parm:0,parmet:9,particular:9,parxial:6,parxialmatrix:5,parxialthinlen:3,pass:[0,5,7,9],path:[0,5],pathelength:[0,5],pathlength:[0,5],pattern:0,paxial:2,peak:9,pencil:[0,4,5],perform:0,pertendicular:6,perzal:2,perzval:2,petzvalsum:2,phase:[0,5,10],phi:8,photop:9,photopicspectrum:1,photot:9,pint:7,pixel:[0,10],place:[2,3,5,6,8,9],plan:[0,4,7],planck:9,planckspectrum:1,plane:[0,2,3,4,5,6],planepair:[2,3],plank:9,plano:2,planoconvex:2,plot:[1,3,4,5,6,9,10],plotpoint:9,plt:3,point:[0,1,2,3,5,6,9],pointinplan:5,poistion:3,polar:8,polynomi:1,pos:[4,5,6],posit:[0,2,3,4,5,6],possibl:2,postit:6,potiton:3,pow:8,power:[3,5,10],prefix:7,present:[5,6,7],princip:[2,3],principl:3,prinicp:3,prinnt:7,print:[1,3,5],printpath:5,priod:9,problem:7,process:[5,7],processd:7,product:8,program:7,progress:5,prompt:[3,7],propag:[3,5,8],propagatethrough:5,propagateto:5,propagationmatrix:3,propat:5,properti:3,propgat:5,provid:7,psf:1,psi:[7,8],pt_or_u:0,pt_or_x:6,pt_or_z:2,pts:4,pupil:2,python:7,quadrat:6,quadric:6,quadricsurfac:1,qudaricsurfac:6,qudric:6,rad:2,radi:6,radial:1,radian:[5,6,7,8],radiu:[0,2,3,4,5,6,10],radu:2,rai:[0,2,4,6,8],random:8,rang:[7,8,9,10],rather:[2,6,7],ratio:[2,6,8],ration:[6,8],raw:10,raymonitor:1,raypath:5,raypenc:0,raypencil:[0,1,4],raytyp:5,read:[1,2,3,9],real:[5,7],recalcul:9,record:[5,7],recoveri:7,rect:8,red:9,redcolourmatch:9,redlimit:9,ref:9,refactiveindex:5,refarct:[3,6,9],refarctiveindex:2,refartiveindex:2,refatc:9,refear:6,refect:8,refer:[0,2,6,9],refix:7,reflect:[5,6,8],reformat:7,refpt:0,refr:[5,6,9],refract:[1,2,3,5,6,8],refractic:8,refractiveindex:[2,5,6,9],refractiveindexcolour:9,refratc:[3,6,9],refratciveindex:5,refrativeindex:[2,6,9],refrctiveindex:2,refrect:9,regular:0,rel:[2,3],releas:7,remad:2,remov:[5,6,7],removeinvalid:5,render:4,repeat:3,repend:2,replac:7,repres:[0,2,3,4,6,9],represent:3,reprompt:7,request:2,requir:[3,7],rerturn:5,reset:[0,2,7,9],resolut:0,respect:3,respons:[7,9],result:[2,9],retail:2,retain:2,retuen:5,retun:8,reurn:2,reus:7,revers:[2,6],rgb:9,right:[3,6],ring:7,root:7,rotat:8,rotateabouti:8,rotateaboutx:8,rotateaboutz:8,rotatian:8,rounc:5,round:[0,8],rubbish:8,s_or_i:5,same:[0,2,3,5,6,9],sampl:0,saniti:7,scale:[2,3,6,8,9],scope:7,scotop:9,scotopicspectrum:1,search:1,second:[2,3,4,8],see:[2,9],select:7,selector:7,self:[0,5,8,9],sensibl:7,separ:3,sepcifi:2,sepectrum:9,serar:3,set:[0,2,3,4,5,6,7,8,9,10],setbend:2,setcentrethick:2,setcurvatur:2,setdefaultwavelength:9,setdegre:8,setedgethick:2,setfixedairindex:9,setfocallength:[2,3],setfromstr:2,setinvalid:[5,8],setiri:2,setjourn:7,setknif:[0,6],setlength:8,setparamet:2,setpoint:[2,6],setpolar:8,setpolardegre:8,setradiu:2,setratio:6,setsiz:6,settemperatur:9,setthick:2,settin:5,setwithpencil:0,setwithrai:[0,4],setwithzernik:0,should:[6,7],side:[2,6],sign:[6,7],signific:3,simpl:[2,7,9],simpler:[2,6],simplesinglet:1,simplest:7,simpli:7,simul:[0,9,10],sinc:[2,3,7],singl:[0,7,9],singlet:1,sinmpl:4,site:9,situ:3,six:2,size:[0,6,10],skew:5,small:[6,7],smooth:9,sourc:[0,2,3,4,5,6,7,8,9,10],sourceplan:5,sourcepoint:[0,1,6],space:[2,3,6,7],special:9,specif:[0,3,4,5,6,7,9],specifi:[0,2,3,4,5,6,7,8,9],specral:9,spectifi:2,spectra:9,spectral:9,spectrum:[1,5],spectum:9,specturm:9,spetrum:5,sphere:8,spheric:6,sphericalsurfac:1,spot:1,spotdiagram:4,spotdigram:1,spread:1,sqr:0,squar:[0,8],standalon:6,standard:[0,2,7,9],start:[0,2,4,5,9],startup:9,state:9,stephen:9,sting:2,stop:7,str:[2,3,4,5,7,9],strike:5,string:[2,3,5,7],structur:7,subset:2,substract:0,success:8,sucess:5,sufac:[6,8],suitabl:5,sum:2,suppi:10,suppl:[8,9],suppli:[0,7],support:[1,3,7,8],sur:5,surafc:6,surafec:2,surfac:[0,1,2,3,5,8],surfaceinteract:[0,6],surfacevector:6,surfc:5,sysout:7,system:[0,2,3],take:[3,6,7],taken:0,taret:0,target:[0,2,3],targetplac:0,targetplan:1,task:0,temperatur:9,tempertaur:9,tempertur:9,tenpertur:9,termin:1,tessar:3,test:[0,5,6,7,8],text:7,than:[2,6,7],thei:[2,3,5],them:[0,3,8],theta:[6,7,8],thi:[0,2,3,4,5,6,7,9,10],thick:[2,3],thickess:2,thicklensmatrix:3,thin:[2,3],thinest:2,thinlensmatrix:3,thinnest:2,though:3,three:[0,2,3,5,6,7,8,9],through:[0,3,5,8],tilt:10,time:5,tio:[1,3],titl:[3,9],tken:2,tprint:7,trace:[0,1,2,3],tragetplan:0,trail:7,transvers:0,trem:0,treturn:3,tri:7,triangl:8,tricolourspectrum:1,trupl:[6,7,8],tupl:[7,8],turn:5,two:[0,2,3,4,5,6,8,9],type:[0,1,2,5,6,8,9,10],typic:[2,3,4,6,7,9],typicv:5,udat:5,unchang:3,under:[0,7],underli:[2,3,6],uniqu:[6,7],unit3d:[0,1,3,5,6],unit:[3,8],uniti:3,unitpair:8,univers:9,unix:7,untest:7,upadt:5,updat:[2,5,9],updatemonitor:5,upper:7,urvatur:2,use:[0,2,5,6,7,8,9],used:[0,2,3,4,5,6,7,8,9],useful:[5,6,7,9],user:[2,3,7,9],usernam:7,uses:7,using:[0,2,3,5,8,9],usual:[2,6,9],util:7,valid:[4,5,6,8,9],valu:[0,1,2,4,5,6,7,8,10],valueerror:6,valul:8,vari:2,variabl:[0,2,5,6,7,9],varial:9,variiabl:9,variou:[0,3,5,6,8],vaue:10,vaul:9,vecftor3d:8,vector2d:[0,1,4,5,6,7,10],vector3:6,vector3d:[0,1,2,3,4,5,6,7],vector:[1,3,5,6],vectors2d:8,vectors3d:8,vectror:8,version:6,vertic:6,vetor3d:8,vextor3d:6,via:[0,3,5,6,7],visual:9,wai:[2,6],wave:[0,2,4,5,9],wave_arrai:9,wavefront:0,wavelelength:4,wavelength:[0,1,2,4,5,6],wavelengthcolour:[5,9],wavelenth:5,wavelnegth:[2,9],wavelngth:2,wavelnngth:2,wavepoint:1,wavepointset:1,wavepooint:0,websit:9,weight:10,what:7,when:[5,7,9],where:[0,2,3,4,5,6,7,8,9],whhere:9,which:[0,2,3,5,6,7,9],whick:2,white:7,whole:[2,5,9],wide:0,width:9,wih:6,within:[8,9],without:[2,6],wold:10,wolf:10,work:[5,7,10],world:7,wrang:9,wright:10,write:7,written:7,wrong:7,wrt:[6,8],wuth:6,www:9,x_or_v:8,xgap:0,xpixel:0,xsize:[0,2,6],xtilt:10,xxx:9,xxxyyi:9,yes:7,ygap:0,you:8,ypixel:0,ysize:[0,2,6],ytilt:10,z_or_v:6,zerkin:10,zern:0,zernik:[0,1],zernikeexpals:1,zernikeexpans:[0,10],zernilkeexpans:0,zero:[0,3,7,8,9],zerom:0,zeromean:0},titles:["Analysis Classes","Welcome to Ray Optics\u2019s documentation!","Lens Classes","Matrix Methods","Point Spread Functions and Spot Diagrams","Ray Tracing","Surface Classes","Terminal Input / Output (tio)","Vector Classes","Wavelength Module","Zernike analysis"],titleterms:{"class":[0,2,3,4,5,6,8,9,10],"default":9,"function":[4,7,9],aberr:0,airindex:9,algebra:3,analysi:[0,10],angl:[7,8],annularapertur:6,basic:7,cauchyindex:9,circularapertur:6,code:3,colour:9,command:7,constant:6,databaselen:2,diagram:4,document:1,doublet:2,escap:7,exampl:[3,7],extend:3,file:7,filenam:7,fixedindex:9,flatsurfac:6,guassianspectrum:9,imageplan:6,index:9,indic:1,infoindex:9,input:7,intens:5,irisapertur:6,journal:7,knifeedgeapertur:6,knifeedgeplot:0,len:2,materialindex:9,matrix:3,menu:7,method:3,modul:9,moment:4,optic:1,opticalgroup:2,opticalimag:0,opticalplan:6,option:7,output:7,parabolicsurfac:6,paraxi:5,paraxialgroup:3,paraxialmatrix:3,paraxialrai:5,photopicspectrum:9,planckspectrum:9,plot:0,point:4,polynomi:10,print:7,psf:4,quadricsurfac:6,radial:10,rai:[1,5],raymonitor:5,raypencil:5,read:7,refract:9,scotopicspectrum:9,simplesinglet:2,singlet:2,sourcepoint:5,spectrum:9,sphericalsurfac:6,spot:4,spotdigram:4,spread:4,support:9,surfac:6,tabl:1,targetplan:0,termin:7,tio:7,trace:5,tricolourspectrum:9,type:7,unit3d:8,valu:9,vector2d:8,vector3d:8,vector:[7,8],wavelength:9,wavepoint:0,wavepointset:0,welcom:1,zernik:10,zernikeexpals:10}})