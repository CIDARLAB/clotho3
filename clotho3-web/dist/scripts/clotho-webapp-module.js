angular.module("clotho.webapp",["clotho.foundation","clotho.interface","clotho.dna","ngSanitize","ngRoute"]),angular.module("clotho.webapp").controller("HomeCtrl",["$scope","Clotho",function(a,b){a.modalContent="<p>Welcome to Clotho!</p><p>Clotho is a platform for automating your genetic engineering projects. Learn how to use Clotho by starting the trail below!</p>",a.enterClotho=function(){b.startTrail("org.clothocad.trails.LearningClotho")},a.enterEugene=function(){b.startTrail("org.clothocad.trails.EugeneCADIntro")},a.enterRaven=function(){b.startTrail("org.clothocad.trails.RavenCADIntro")}}]),angular.module("clotho.webapp").controller("TeamCtrl",["$scope","Clotho",function(a,b){b.query({id:{$regex:"clotho.developer.*",$options:"i"}},{mute:!0}).then(function(b){a.ClothoTeam=b})}]),angular.module("clotho.webapp").controller("BrowserCtrl",["$scope","Clotho","$filter",function(a,b,c){b.recent().then(function(b){a.recent_array=b,a.sort(!1)}),a.nav={user:{text:"My Stuff",subtext:"Your Digital Locker",value:"user"},group:{text:"Our Group",subtext:"Your Group's Info",value:"group"},all:{text:"Everything",subtext:"All of Clotho",value:"all"}},a.setCurrent=function(b){a.current=b},a.sort=function(b){if(b){if(a.catSort)return;a.catSort=!0,a.recent=c("categorize")(a.recent_array,"type"),a.recent.Instance=c("categorize")(a.recent.Instance,"schema.name")}else{if(!a.catSort&&angular.isDefined(a.catSort))return;a.catSort=!1,a.recent={entries:angular.copy(a.recent_array)}}},a.base64icon=base64icon}]);var base64icon="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEAAAABACAYAAACqaXHeAAACtElEQVR4Xu2Y3UtqURDFlyYVqQhBQon5koqYiUIElSD951qa+AkGUWDUgxC9KEhqfnfXgOLl3gKP2Xlw5kXO0T2zZ+2Z+eG2NBqNCdbYLCqAVoC2gM6ANZ6B0CGoFFAKKAWUAkqBNVZAMagYVAwqBhWDawwB/TOkGFQMKgYVg4pBxeAaK7A0Bh8fH/H6+orJZIL9/X0Eg0FYLJaZpH8og/v7e+zs7CAWi/313Ve6r8LnV7GWEuDh4QH1eh02m038D4dD+Hw++P1+eaYoqVRK3m9ubuLy8hJWq/XbeluFz+8CGhZgNBohnU7LiSYSCQwGA7y8vMDpdMLj8UjMWq0m72jb29u4uLhAu91GtVrFxsYG4vE43t/fwaS3trYQiURwc3OzsM/5ilu0mw0L0O/3ZbPj8ViS/vj4kBYIBAKyh263i0wmg8PDQ7RaLXQ6HRGKmy2VSmg2m9jd3ZV1/I5Vw/VGfS6a+PT3hgVgb5fLZfHDFmCZ05hEOByWJJl4MplEPp9Hr9ebCTAvHte4XC6cnp5iGZ+/LgBPOJvNSumen5/LKeZyOen1k5MTFItF2O12eL1ePD09SaUcHR3JM+35+Vne087OzqSKlvVpRATDFcAT5wxgb1OA6TMFCIVCqFQq/+yHA/Dq6gqcH9fX17Oq2dvbQzQanfkw4tNI8lxjWABOeFYAT+3g4EA+2ddutxvHx8dgmdOYdKFQkOR40kzu7u4Ob29vcDgcUjmsDg5ArjXq89cFYEAOsNvbWzlRGtuBvcwk542tQQFIAc4FCjKlB4Ug9zlHpjRZ1KcpFJhPkCij8UR/ylbh8397M9wCP5Wo2X5UAL0R0hshvRHSGyGzJ7GZ8ZUCSgGlgFJAKWDmFDY7tlJAKaAUUAooBcyexGbGVwooBZQCSgGlgJlT2OzYSoF1p8AnDSiNnx2jBucAAAAASUVORK5CYII=";angular.module("clotho.webapp").controller("EditorCtrl",["$scope","$route","$location","Clotho","ClothoSchemas",function(a,b,c,d,e){a.$watch("editable.id",function(a){c.search("id",a||null)}),a.$on("$routeUpdate",function(b,c){var d=c.params.id;angular.isEmpty(a.editable)||a.editable.id==d||(a.editableId=d)}),a.$on("$destroy",function(){c.search("id",null).replace()}),b.current.params.id&&(a.editableId=b.current.params.id),a.editModePass=!1,a.objectTypes=e.sharableTypes,a.schemas=[],e.retrievedSchemas.then(function(b){a.schemas=b}),a.queryWrapper=function(a){return d.autocomplete(a).then(function(a){return a||[]})},a.createNewNonInstance=function(b){a.editable=e.createScaffold(b),a.editModePass=!0},a.createNewInstance=function(b,c){a.editable=e.createScaffold(c),a.editModePass=!0},a.editExisting=function(b,c){d.get(c).then(function(b){a.editable=b,a.editModePass=!0},function(){console.log("for some reason could not get that sharable...")})}}]),angular.module("clotho.trails").controller("TrailCtrl",["$scope","$route","$timeout","Clotho","Trails","hotkeys","$location",function(a,b,c,d,e,f,g){a.id=b.current.params.id,a.trail=b.current.locals.trail,b.current.params.position&&(a.current=b.current.params.position,a.currentPage=e.extractPage(a.trail,a.current)),a.activate=e.activate,a.favorite=function(){e.favorite(a.id)},a.share=function(){e.share(a.id)},a.next=function(){a.activate(e.calcNextPage(a.trail,a.current))},a.prev=function(){a.activate(e.calcPrevPage(a.trail,a.current))},a.mapIcon=e.mapIcon,a.$on("$routeUpdate",function(b,c){a.current!=c.params.position&&(a.current=c.params.position||null,a.currentPage=e.extractPage(a.trail,a.current),angular.isEmpty(a.currentPage)&&a.activate("0-0"))}),a.$on("$destroy",function(){g.search("position",null).replace(),g.search("id",null).replace()}),a.activate(a.current)}]),angular.module("clotho.webapp").controller("TrailsCtrl",["$scope","$location","Clotho",function(a,b,c){a.headerMockTrail={name:"Trail Browser",description:"Learn Synthetic Biology with Clotho"},a.trails=[],c.query({schema:"org.clothocad.model.Trail"}).then(function(b){a.trails=b}),a.defaultTrailIcon="images/trails/trails_logo.png",a.startTrail=function(a){c.startTrail(a)},a.startTrailPage=function(a,d){b.search("position",d),c.startTrail(a)},a.highlight=function(b){a.highlighted=b,a.loading=!0,c.get(b.id,{mute:!0}).then(function(b){a.loading=!1,a.selected=b})}}]),angular.module("clotho.trails").controller("TrailSplashCtrl",["$scope","$location","Clotho",function(a,b,c){a.topics=[{title:"Introducing Clotho",trails:[{title:"Introduction to Clotho",id:"org.clothocad.trails.LearningClotho"}]},{title:"Basic Molecular Biology",trails:[{title:"Introduction to Synthetic Biology",id:"org.clothocad.trails.youtube.IntroductiontoSyntheticBiology"},{title:"Organic Chemistry and Molecular Biology",id:"org.clothocad.trails.youtube.OrganicChemistryandMolecularBiology"},{title:"Basics of OOP and Organic Chemistry",id:"org.clothocad.trails.youtube.BasicsofOOPandOrganicChemistry"}]},{title:"Fabrication",trails:[{title:"DNA Manipulation Enzymes",id:"org.clothocad.trails.youtube.DNAManipulationEnzymes"},{title:"DNA Fabrication",id:"org.clothocad.trails.youtube.DNAFabrication"},{title:"Genome Manipulation",id:"org.clothocad.trails.youtube.GenomeManipulation"},{title:"Chassis and Strains",id:"org.clothocad.trails.youtube.ChassisandStrains"},{title:"Parts",id:"org.clothocad.trails.youtube.Parts"},{title:"Human Practices",id:"org.clothocad.trails.youtube.HumanPractices"}]},{title:"Analysis",trails:[{title:"Whole-Cell Analysis",id:"org.clothocad.trails.youtube.Whole-CellAnalysis"},{title:"In vitro Analysis",id:"org.clothocad.trails.youtube.InvitroAnalysis"},{title:"Directed Evolution",id:"org.clothocad.trails.youtube.DirectedEvolution"},{title:"Combinatorial Libraries",id:"org.clothocad.trails.youtube.CombinatorialLibraries"},{title:"Automation Tools",id:"org.clothocad.trails.youtube.AutomationTools"}]},{title:"Design",trails:[{title:"Introduction to Biosynthesis",id:"org.clothocad.trails.youtube.IntroductiontoBiosynthesis"},{title:"Literature Search and Pathway Discovery",id:"org.clothocad.trails.youtube.LiteratureSearchandPathwayDiscovery"},{title:"Metabolic Optimization",id:"org.clothocad.trails.youtube.MetabolicOptimization"},{title:"Monofunctional Enzyme Biosynthesis Examples ",id:"org.clothocad.trails.youtube.MonofunctionalEnzymeBiosynthesisExamples"},{title:"Polymeric Metabolites",id:"org.clothocad.trails.youtube.PolymericMetabolites"},{title:"Modular PKS and NRP Engineering",id:"org.clothocad.trails.youtube.ModularPKSandNRPEngineering"},{title:"Processes and Localization",id:"org.clothocad.trails.youtube.ProcessesandLocalization"},{title:"Transcription Models",id:"org.clothocad.trails.youtube.TranscriptionModels"},{title:"Translation Models",id:"org.clothocad.trails.youtube.TranslationModels"},{title:"Biochemical Circuits",id:"org.clothocad.trails.youtube.BiochemicalCircuits"},{title:"Stochastic modeling",id:"org.clothocad.trails.youtube.Stochasticmodeling"},{title:"Biological Relationships",id:"org.clothocad.trails.youtube.BiologicalRelationships"},{title:"Engineering microbial interactions",id:"org.clothocad.trails.youtube.Engineeringmicrobialinteractions"},{title:"Refactoring",id:"org.clothocad.trails.youtube.Refactoring"},{title:"Advanced microbial processes",id:"org.clothocad.trails.youtube.Advancedmicrobialprocesses"},{title:"Eukaryotic Devices",id:"org.clothocad.trails.youtube.EukaryoticDevices"},{title:"Developmental Devices",id:"org.clothocad.trails.youtube.DevelopmentalDevices"},{title:"Part Engineering",id:"org.clothocad.trails.youtube.PartEngineering"}]},{title:"Wetlab: Cloning Experiment",trails:[{title:"Overview of Basic Cloning",id:"org.clothocad.trails.youtube.OverviewofBasicCloning"},{title:"Running PCR",id:"org.clothocad.trails.youtube.RunningPCR"},{title:"Digesting DNAs with restriction endonucleases",id:"org.clothocad.trails.youtube.DigestingDNAswithrestrictionendonucleases"},{title:"Gel Purification",id:"org.clothocad.trails.youtube.GelPurification"},{title:"Ligation and Transformation",id:"org.clothocad.trails.youtube.LigationandTransformation"},{title:"Growing up colonies and then prepping DNA",id:"org.clothocad.trails.youtube.GrowingupcoloniesandthenpreppingDNA"}]},{title:"Advanced Tools",trails:[]}],a.startTrail=function(a){console.log(a),c.startTrail(a)},a.startTrailPage=function(a,d){b.search("position",d),c.startTrail(a)},a.highlight=function(b){a.highlighted=b,a.loading=!0,c.get(b.id,{mute:!0}).then(function(b){a.loading=!1,a.selected=b})}}]),angular.module("clotho.webapp").controller("WidgetsCtrl",["$scope","$compile",function(a,b){a.myObj={myProp:"aacccggttt"},a.myModel="WHATS up",a.myDNAmodel="aaaagt",a.bootstrapCallback=function(a){console.log("widget controller callback! passed element:",a)},a.bootstrapNewApp=function(){console.log("bootstrap start");var c=angular.element('<clotho-show id="123456789" callback="bootstrapCallback"></clotho-show>');angular.element(document).find("insert-widgets-here").append(c),b(c)(a),console.log("done compiling")},a.bootstrapNewAppExternal=function(){var c=angular.element('<clotho-show id="123456789"></clotho-show>');angular.element(document).find("#clothoAppWidgets").append(c),b(c)(a)},a.bootstrapSimple=function(){var c=angular.element('<clotho-show id="813579135"></clotho-show>');angular.element(document).find("insert-widgets-here").append(c),b(c)(a)},a.someValue="from parent controller not passed down",a.openCallback=function(){console.log("modal openeed")},a.closeCallback=function(){console.log("modal closed")}}]);