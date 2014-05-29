"use strict";angular.module("clotho.webapp").controller("TestSchemaviewCtrl",["$scope","Clotho","$clothoModal",function(a,b,c){b.query({name:"maxbates"}).then(function(b){a.retrievedId=b.length?b[0].id:""}),a.createModal=function(){var b=a.$new();b.modalActions=[{"class":"info",text:"Great!",action:c.destroy}],c.create({title:"Hey there",content:"'here is some <br>content'",actions:"modalActions"},b)}}]),angular.module("clotho.webapp").controller("TestPlaylistimportCtrl",["$scope","Youtube","Clotho","$q",function(a,b){a.playlistId="PL2aPXzks-TgO0k9PhT__NSh2x6HNimaOy",a.$watch("playlistId",function(c){b.playlistItems(c).then(function(b){a.playlistInfo=b}),b.playlistInfo(c).then(function(b){a.playlistItems=b}),b.playlistToTrail(c).then(function(b){a.playlistTrail=b})}),a.$watch("search",function(c){c?b.videoSearch(c).then(function(b){a.searchResult=b}):a.searchResult=""})}]),angular.module("clotho.webapp").controller("TestQuizCtrl",["$scope","$timeout",function(a,b){b(function(){a.trueFalse={question:{type:"truefalse",title:"True False",question:"Clotho is great?"},grade:{answer:{type:"boolean",value:!0}}}},1e3),a.number={question:{type:"number",title:"Number with tolerance (0.01)",question:"What is Pi?"},grade:{answer:{type:"number",tolerance:.01,value:3.1415926}}},a.staticMC={question:{type:"mc",title:"Static Multiple Choice",question:"Find the reverse complement of the sequence ACCGGGTTTT",options:["accgggtttt","TGGCCCAAAA","TTTTGGGCCA","AAAACCCGGT"],hint:"Example hint!"},options:{showAnswer:!0,allowMultiple:!0,allowRetry:!0,randomization:!1},grade:{answer:{type:"string",value:"AAAACCCGGT"}}},a.staticTemplating={question:{type:"mc",title:"Static Templating",question:"Find the reverse complement of the sequence {{mySeq}}",options:["{{value1}}","{{value2}}","{{value3}}","{{value4}}"]},dictionary:{"static":{mySeq:"ACCGGGTTTT",value1:"accgggtttt",value2:"TGGCCCAAAA",value3:"TTTTGGGCCA",value4:"AAAACCCGGT"}},grade:{answer:{type:"string",value:"AAAACCCGGT"}}},a.staticFeedback={question:{type:"mc",title:"Static Feedback",question:"Find the reverse complement of the sequence ACCGGGTTTT",options:["accgggtttt","TGGCCCAAAA","TTTTGGGCCA","AAAACCCGGT"]},options:{allowRetry:!0},grade:{answer:{type:"string",value:"AAAACCCGGT"}},feedback:{"default":"To reverse complement a sequence, simply reverse the order, and complement each base (order in which you do this does not matter)","static":{accgggtttt:"Nope, that's the lowercase",TGGCCCAAAA:"Nope, that's the complement",TTTTGGGCCA:"Nope, that's the reverse"}}},a.staticFillin={question:{type:"fillin",title:"Static Fillin",question:"Find the reverse complement of the sequence AAACCGT"},options:{showAnswer:!0},grade:{answer:{type:"string",value:"ACGGTTT"}}},a.dynamicTemplating={question:{type:"mc",title:"Dynamic Templating",question:"Find the reverse complement of the sequence {{mySeq}}",options:["{{value1}}","{{value2}}","{{value3}}","{{mySeq}}"],hint:"When you retry, all templated values are recalculated. They are then passed to the server for grading"},options:{allowRetry:!0},dictionary:{dynamic:[{mySeq:{id:"org.clothocad.test.randomSequence",args:[16]}},{value1:{id:"org.clothocad.test.revcomp",args:["{{mySeq}}"]}},{value2:{id:"org.clothocad.test.reverse",args:["{{mySeq}}"]}},{value3:{id:"org.clothocad.test.complement",args:["{{mySeq}}"]}}]},grade:{answer:{type:"function",value:"org.clothocad.test.revcomp"},args:["{{mySeq}}"]}},a.functionGrading={question:{type:"fillin",title:"Own grading function",question:"Copy in this sequence: {{mySeq}}",hint:"This question has its own grading function, for more complicated logic. Arguments can be specified by author. This one just makes sure args + input are equal. This function is not case sensitive (but grading normally is)"},options:{allowRetry:!0},dictionary:{dynamic:[{mySeq:{id:"org.clothocad.test.randomSequence",args:[40]}}]},grade:{"function":"org.clothocad.test.customGradeFunction",args:["{{mySeq}}"]}}}]),angular.module("clotho.webapp").controller("TestConstructionCtrl",["$scope","$http","ConstructionSimulator",function(a,b,c){b.get("models/construction/construction_pHA581.json").success(function(b){a.demoConstruction=b,c.process(b).then(function(b){a.processed=b})}),b.get("models/construction/construction_parsed_kan.json").success(function(b){a.parsed=b})}]),angular.module("clotho.trails").controller("TestContstructiontrailCtrl",["$scope","$route","$timeout","Clotho","Trails","$location",function(a,b,c,d,e,f){a.id=b.current.params.id,a.trail=b.current.locals.trail,b.current.params.position&&(a.current=b.current.params.position,a.currentPage=e.extractPage(a.trail,a.current)),a.activate=e.activate,a.favorite=function(){e.favorite(a.id)},a.share=function(){e.share(a.id)},a.next=function(){a.activate(e.calcNextPage(a.trail,a.current))},a.prev=function(){a.activate(e.calcPrevPage(a.trail,a.current))},a.mapIcon=e.mapIcon,a.$on("$routeUpdate",function(b,c){a.current!=c.params.position&&(a.current=c.params.position||null,a.currentPage=e.extractPage(a.trail,a.current),angular.isEmpty(a.currentPage)&&a.activate("0-0"))}),a.$on("$destroy",function(){f.search("position",null).replace(),f.search("id",null).replace()}),a.activate(a.current)}]),angular.module("clotho.webapp").controller("TestFocusCtrl",["$scope","CommandBar","$focus","$parse",function(a,b,c){a.setInput=function(){b.setInput("revcomp aacc")},a.type=function(){c.typeOutSearch("revcomp aacc")}}]);