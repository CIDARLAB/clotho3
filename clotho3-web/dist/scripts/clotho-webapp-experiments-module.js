"use strict";angular.module("clotho.webapp").controller("TestTokenizerCtrl",["$scope",function(a){a.myModel=[]}]),angular.module("clotho.webapp").controller("TestSchemaviewCtrl",["$scope","Clotho","$clothoModal",function(a,b,c){b.query({name:"maxbates"}).then(function(b){a.retrievedId=b.length?b[0].id:""}),a.createModal=function(){var b=a.$new();b.modalActions=[{"class":"info",text:"Great!",action:c.destroy}],c.create({title:"Hey there",content:"'here is some <br>content'",actions:"modalActions"},b)}}]),angular.module("clotho.trails").controller("TestTrailCtrl",["$scope","$route","$timeout","Clotho","Trails","hotkeys","$location",function(a,b,c,d,e,f,g){a.id=b.current.params.id,a.trail=b.current.locals.trail,b.current.params.position&&(a.current=b.current.params.position,a.currentPage=e.extractPage(a.trail,a.current)),a.activate=e.activate,a.favorite=function(){e.favorite(a.id)},a.share=function(){e.share(a.id)},a.next=function(){a.activate(e.calcNextPage(a.trail,a.current))},a.prev=function(){a.activate(e.calcPrevPage(a.trail,a.current))},a.mapIcon=e.mapIcon,a.$on("$routeUpdate",function(b,c){a.current!=c.params.position&&(a.current=c.params.position||null,a.currentPage=e.extractPage(a.trail,a.current),angular.isEmpty(a.currentPage)&&a.activate("0-0"))}),a.$on("$destroy",function(){g.search("position",null).replace(),g.search("id",null).replace()}),a.activate(a.current)}]),angular.module("clotho.trails").controller("TestTrailBrowserCtrl",["$scope","Clotho",function(a,b){b.query({schema:"org.clothocad.model.Trail"}).then(function(b){a.trails=b}),a.startTrail=b.startTrail}]),angular.module("clotho.webapp").controller("TestPlaylistimportCtrl",["$scope","Youtube","Clotho","$q",function(a,b){a.playlistId="PL2aPXzks-TgO0k9PhT__NSh2x6HNimaOy",a.$watch("playlistId",function(c){b.playlistItems(c).then(function(b){a.playlistInfo=b}),b.playlistInfo(c).then(function(b){a.playlistItems=b}),b.playlistToTrail(c).then(function(b){a.playlistTrail=b})}),a.$watch("search",function(c){c?b.videoSearch(c).then(function(b){a.searchResult=b}):a.searchResult=""})}]),angular.module("clotho.webapp").controller("TestQuizCtrl",["$scope","$timeout",function(a,b){b(function(){a.trueFalse={question:{type:"truefalse",title:"True False",question:"Clotho is great?"},grade:{answer:{type:"boolean",value:!0}}}},1e3),a.number={question:{type:"number",title:"Number with tolerance (0.01)",question:"What is Pi?"},grade:{answer:{type:"number",tolerance:.01,value:3.1415926}}},a.staticMC={question:{type:"mc",title:"Static Multiple Choice",question:"Find the reverse complement of the sequence ACCGGGTTTT",options:["accgggtttt","TGGCCCAAAA","TTTTGGGCCA","AAAACCCGGT"],hint:"Example hint!"},options:{showAnswer:!0,allowMultiple:!0,allowRetry:!0,randomization:!1},grade:{answer:{type:"string",value:"AAAACCCGGT"}}},a.staticTemplating={question:{type:"mc",title:"Static Templating",question:"Find the reverse complement of the sequence {{mySeq}}",options:["{{value1}}","{{value2}}","{{value3}}","{{value4}}"]},dictionary:{"static":{mySeq:"ACCGGGTTTT",value1:"accgggtttt",value2:"TGGCCCAAAA",value3:"TTTTGGGCCA",value4:"AAAACCCGGT"}},grade:{answer:{type:"string",value:"AAAACCCGGT"}}},a.staticFeedback={question:{type:"mc",title:"Static Feedback",question:"Find the reverse complement of the sequence ACCGGGTTTT",options:["accgggtttt","TGGCCCAAAA","TTTTGGGCCA","AAAACCCGGT"]},options:{allowRetry:!0},grade:{answer:{type:"string",value:"AAAACCCGGT"}},feedback:{"default":"To reverse complement a sequence, simply reverse the order, and complement each base (order in which you do this does not matter)","static":{accgggtttt:"Nope, that's the lowercase",TGGCCCAAAA:"Nope, that's the complement",TTTTGGGCCA:"Nope, that's the reverse"}}},a.staticFillin={question:{type:"fillin",title:"Static Fillin",question:"Find the reverse complement of the sequence AAACCGT"},options:{showAnswer:!0},grade:{answer:{type:"string",value:"ACGGTTT"}}},a.dynamicTemplating={question:{type:"mc",title:"Dynamic Templating",question:"Find the reverse complement of the sequence {{mySeq}}",options:["{{value1}}","{{value2}}","{{value3}}","{{mySeq}}"],hint:"When you retry, all templated values are recalculated. They are then passed to the server for grading"},options:{allowRetry:!0},dictionary:{dynamic:[{mySeq:{id:"org.clothocad.test.randomSequence",args:[16]}},{value1:{id:"org.clothocad.test.revcomp",args:["{{mySeq}}"]}},{value2:{id:"org.clothocad.test.reverse",args:["{{mySeq}}"]}},{value3:{id:"org.clothocad.test.complement",args:["{{mySeq}}"]}}]},grade:{answer:{type:"function",value:"org.clothocad.test.revcomp"},args:["{{mySeq}}"]}},a.functionGrading={question:{type:"fillin",title:"Own grading function",question:"Copy in this sequence: {{mySeq}}",hint:"This question has its own grading function, for more complicated logic. Arguments can be specified by author. This one just makes sure args + input are equal. This function is not case sensitive (but grading normally is)"},options:{allowRetry:!0},dictionary:{dynamic:[{mySeq:{id:"org.clothocad.test.randomSequence",args:[40]}}]},grade:{"function":"org.clothocad.test.customGradeFunction",args:["{{mySeq}}"]}}}]),angular.module("clotho.trails").controller("TestTrailSplashCtrl",["$scope","Clotho","Youtube",function(a,b){a.topics=[{title:"Introducing Clotho",trails:[{title:"Introduction to Clotho",id:"org.clothocad.trails.LearningClotho"}]},{title:"Basic Molecular Biology",trails:[{title:"Introduction to Synthetic Biology",id:"org.clothocad.trails.youtube.IntroductiontoSyntheticBiology"},{title:"Organic Chemistry and Molecular Biology",id:"org.clothocad.trails.youtube.OrganicChemistryandMolecularBiology"}]},{title:"Fabrication",trails:[{title:"DNA Manipulation Enzymes",id:"org.clothocad.trails.youtube.DNAManipulationEnzymes"},{title:"DNA Fabrication",id:"org.clothocad.trails.youtube.DNAFabrication"},{title:"Genome Manipulation",id:"org.clothocad.trails.youtube.GenomeManipulation"},{title:"Chassis and Strains",id:"org.clothocad.trails.youtube.ChassisandStrains"},{title:"Parts",id:"org.clothocad.trails.youtube.Parts"},{title:"Human Practices",id:"org.clothocad.trails.youtube.HumanPractices"}]},{title:"Analysis",trails:[{title:"Whole-Cell Analysis",id:"org.clothocad.trails.youtube.Whole-CellAnalysis"},{title:"In vitro Analysis",id:"org.clothocad.trails.youtube.InvitroAnalysis"},{title:"Directed Evolution",id:"org.clothocad.trails.youtube.DirectedEvolution"},{title:"Combinatorial Libraries",id:"org.clothocad.trails.youtube.CombinatorialLibraries"},{title:"Automation Tools",id:"org.clothocad.trails.youtube.AutomationTools"}]},{title:"Design",trails:[{title:"Introduction to Biosynthesis",id:"org.clothocad.trails.youtube.IntroductiontoBiosynthesis"},{title:"Literature Search and Pathway Discovery",id:"org.clothocad.trails.youtube.LiteratureSearchandPathwayDiscovery"},{title:"Metabolic Optimization",id:"org.clothocad.trails.youtube.MetabolicOptimization"},{title:"Monofunctional Enzyme Biosynthesis Examples ",id:"org.clothocad.trails.youtube.MonofunctionalEnzymeBiosynthesisExamples"},{title:"Polymeric Metabolites",id:"org.clothocad.trails.youtube.PolymericMetabolites"},{title:"Modular PKS and NRP Engineering",id:"org.clothocad.trails.youtube.ModularPKSandNRPEngineering"},{title:"Processes and Localization",id:"org.clothocad.trails.youtube.ProcessesandLocalization"},{title:"Transcription Models",id:"org.clothocad.trails.youtube.TranscriptionModels"},{title:"Translation Models",id:"org.clothocad.trails.youtube.TranslationModels"},{title:"Biochemical Circuits",id:"org.clothocad.trails.youtube.BiochemicalCircuits"},{title:"Stochastic modeling",id:"org.clothocad.trails.youtube.Stochasticmodeling"},{title:"Biological Relationships",id:"org.clothocad.trails.youtube.BiologicalRelationships"},{title:"Engineering microbial interactions",id:"org.clothocad.trails.youtube.Engineeringmicrobialinteractions"},{title:"Refactoring",id:"org.clothocad.trails.youtube.Refactoring"},{title:"Advanced microbial processes",id:"org.clothocad.trails.youtube.Advancedmicrobialprocesses"},{title:"Eukaryotic Devices",id:"org.clothocad.trails.youtube.EukaryoticDevices"},{title:"Developmental Devices",id:"org.clothocad.trails.youtube.DevelopmentalDevices"},{title:"Part Engineering",id:"org.clothocad.trails.youtube.PartEngineering"}]},{title:"Wetlab: Cloning Experiment",trails:[{title:"Overview of Basic Cloning",id:"org.clothocad.trails.youtube.OverviewofBasicCloning"},{title:"Running PCR",id:"org.clothocad.trails.youtube.RunningPCR"},{title:"Digesting DNAs with restriction endonucleases",id:"org.clothocad.trails.youtube.DigestingDNAswithrestrictionendonucleases"},{title:"Gel Purification",id:"org.clothocad.trails.youtube.GelPurification"},{title:"Ligation and Transformation",id:"org.clothocad.trails.youtube.LigationandTransformation"},{title:"Growing up colonies and then prepping DNA",id:"org.clothocad.trails.youtube.GrowingupcoloniesandthenpreppingDNA"}]},{title:"Advanced Tools",trails:[]}],a.startTrail=b.startTrail,a.highlight=function(c){a.highlighted=c,a.loading=!0,b.get(c.id).then(function(b){a.loading=!1,a.selected=b})}}]),angular.module("clotho.webapp").controller("TestConstructionCtrl",["$scope","$http",function(a,b){b.get("models/construction/construction_parsed_kan.json").success(function(b){a.parsed=b})}]),angular.module("clotho.webapp").controller("TestContstructiontrailCtrl",["$scope","$route","$timeout","Clotho","Trails","$location",function(a,b,c,d,e,f){a.id=b.current.params.id,a.trail=b.current.locals.trail,b.current.params.position&&(a.current=b.current.params.position,a.currentPage=e.extractPage(a.trail,a.current)),a.activate=e.activate,a.favorite=function(){e.favorite(a.id)},a.share=function(){e.share(a.id)},a.next=function(){a.activate(e.calcNextPage(a.trail,a.current))},a.prev=function(){a.activate(e.calcPrevPage(a.trail,a.current))},a.mapIcon=e.mapIcon,a.$on("$routeUpdate",function(b,c){a.current!=c.params.position&&(a.current=c.params.position||null,a.currentPage=e.extractPage(a.trail,a.current),angular.isEmpty(a.currentPage)&&a.activate("0-0"))}),a.$on("$destroy",function(){f.search("position",null).replace(),f.search("id",null).replace()}),a.activate(a.current)}]);