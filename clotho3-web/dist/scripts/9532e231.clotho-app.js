"use strict";angular.module("clotho.fullPackage",["clotho.foundation","clotho.commandbar","clotho.webapp","clotho.editor","clotho.interface","clotho.trails","clotho.dna","clotho.construction"]),angular.module("clothoRoot",["clotho.fullPackage"]).config(["$routeProvider","$locationProvider",function(a,b){b.html5Mode(!1).hashPrefix("!"),a.when("/",{templateUrl:"views/home.html",controller:"HomeCtrl",title:"Home",hotkeys:[["h","Show Intro Modal","showHelp = !showHelp"]]}).when("/settings",{templateUrl:"views/settings.html",controller:"SettingsCtrl",title:"Settings"}).when("/about",{templateUrl:"views/about.html",title:"About"}).when("/team",{templateUrl:"views/team.html",controller:"TeamCtrl",title:"Team"}).when("/browser",{templateUrl:"views/browser.html",controller:"BrowserCtrl",title:"Browser"}).when("/edit",{redirectTo:"/editor"}).when("/editor",{templateUrl:"views/editor.html",controller:"EditorCtrl",title:"Editor",reloadOnSearch:!1,resolve:{deps:["codemirrorLoader",function(a){return a.loadMain()}]}}).when("/executor",{templateUrl:"views/executor.html",controller:"ExecutorCtrl",title:"Function Executor",reloadOnSearch:!1}).when("/trails",{templateUrl:"views/trails.html",controller:"TrailsCtrl",title:"Trails"}).when("/trail",{templateUrl:"views/trail.html",controller:"TrailCtrl",title:"Trail",reloadOnSearch:!1,resolve:{trail:["Clotho","$q","$http","$route","$location","Trails",function(a,b,c,d,e,f){if(angular.isUndefined(d.current.params.id))return b.when(null);var g=b.defer();return a.get(d.current.params.id).then(function(a){f.compile(a).then(function(b){d.current.$$route.title=a.name,g.resolve(b)})},function(){e.path("/trails")}),g.promise}]},hotkeys:[["alt+left","Previous page of Trail","prev()"],["alt+right","Next page of Trail","next()"]]}).when("/trail-splash",{templateUrl:"views/trail-splash.html",controller:"TrailSplashCtrl",title:"Trail Splash"}).when("/terminal",{templateUrl:"views/_command/terminal.html",title:"Terminal",resolve:{deps:function(){return $clotho.extensions.mixin("scripts/_command/terminal.js")}}}).when("/import",{templateUrl:"views/import/intro.html",controller:"ImportCtrl",title:"Import Wizard"}).when("/import/youtubePlaylist",{templateUrl:"views/import/youtubePlaylist.html",controller:"YoutubePlaylistImportCtrl"}).when("/import/ape",{templateUrl:"views/import/ape.html",controller:"ApeImportCtrl",title:"GenBank Import"}).when("/import/ncbi",{templateUrl:"views/import/ncbi.html",controller:"NCBIImportCtrl",title:"NCBI Import"}).when("/widgets",{templateUrl:"views/widgets.html",controller:"WidgetsCtrl"}).when("/test/schemaview",{templateUrl:"views/test/schemaview.html",controller:"TestSchemaviewCtrl",hotkeys:[["m","Show Programmatic Modal","createModal()"]]}).when("/test/quiz",{templateUrl:"views/test/quiz.html",controller:"TestQuizCtrl"}).when("/test/construction",{templateUrl:"views/test/construction.html",controller:"TestConstructionCtrl",title:"Construction Test"}).when("/test/constructionTrail",{templateUrl:"views/test/contstructiontrail.html",controller:"TestContstructiontrailCtrl",title:"Construction Trail Test",reloadOnSearch:!1,resolve:{trail:["$q","$http","$route","Trails",function(a,b,c,d){var e=a.defer();return b.get("models/org.clothocad.trails.constructionFiles.json").then(function(a){d.compile(a.data).then(function(a){e.resolve(a)})}),e.promise}]}}).when("/test/focus",{templateUrl:"views/test/focus.html",controller:"TestFocusCtrl"}).otherwise({redirectTo:"/"})}]).run(["$rootScope","$route","$window","interfaceConfig",function(a,b,c,d){a.$on("$routeChangeSuccess",function(a,b){var d=angular.isDefined(b.$$route)?b.$$route.title:null;c.document.title="Clotho"+(angular.isDefined(d)?" | "+d:"")}),a.interfaceConfig=d}]);