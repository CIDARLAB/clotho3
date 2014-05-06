"use strict";angular.module("clotho.fullPackage",["clotho.foundation","clotho.commandbar","clotho.webapp","clotho.editor","clotho.interface","clotho.trails"]),angular.module("clothoRoot",["clotho.fullPackage"]).config(["$routeProvider",function(a){a.when("/",{templateUrl:"views/home.html",controller:"HomeCtrl",title:"Home",hotkeys:[["h","Show Intro Modal","showHelp = !showHelp"]]}).when("/about",{templateUrl:"views/about.html",controller:"AboutCtrl"}).when("/team",{templateUrl:"views/team.html",controller:"TeamCtrl"}).when("/browser",{templateUrl:"views/browser.html",controller:"BrowserCtrl"}).when("/edit",{redirectTo:"/editor"}).when("/editor",{templateUrl:"views/editor.html",controller:"EditorCtrl",reloadOnSearch:!1,resolve:{deps:["codemirrorLoader",function(a){return a.loadMain()}]}}).when("/trails",{templateUrl:"views/trails.html",controller:"TrailsCtrl"}).when("/trails/:id",{templateUrl:"views/trail.html",controller:"TrailCtrl",reloadOnSearch:!1,resolve:{trail:["Clotho","$q","$http","$route","Trails",function(a,b,c,d,e){var f=b.defer();return a.get(d.current.params.id).then(function(a){e.compile(a).then(function(b){d.current.$$route.title=a.name,f.resolve(b)})}),f.promise}]},hotkeys:[["alt+left","Previous page of Trail","prev()"],["alt+right","Next page of Trail","next()"]]}).when("/terminal",{templateUrl:"views/_command/terminal.html",title:"Terminal",resolve:{deps:function(){return $clotho.extensions.mixin("scripts/_command/terminal.js")}}}).when("/widgets",{templateUrl:"views/widgets.html",controller:"WidgetsCtrl"}).when("/test/tokenizer",{templateUrl:"views/test/tokenizer.html",controller:"TestTokenizerCtrl"}).when("/test/schemaview",{templateUrl:"views/test/schemaview.html",controller:"TestSchemaviewCtrl",hotkeys:[["m","Show Programmatic Modal","createModal()"]]}).when("/test/trail",{templateUrl:"views/test/trail.html",controller:"TestTrailCtrl"}).when("/test/trail-browser",{templateUrl:"views/test/trail-browser.html",controller:"TestTrailBrowserCtrl"}).when("/test/trail-overview",{templateUrl:"views/test/trail-overview.html",controller:"TestTrailOverviewCtrl"}).when("/test/playlistimport",{templateUrl:"views/test/playlistimport.html",controller:"TestPlaylistimportCtrl"}).otherwise({redirectTo:"/"})}]).run(["$rootScope","$route","$window","$location","Clotho","CommandBar","hotkeys",function(a,b,c,d,e,f,g){a.$on("$routeChangeSuccess",function(a,b){var d=b.$$route.title;c.document.title="Clotho"+(angular.isDefined(d)?" | "+d:"")}),g.add("f","Focus Command Bar",function(a){a.preventDefault(),f.focusInput()}),g.add("a","Show Activity Log",function(a){a.preventDefault(),f.showActivityLog()}),g.add("g h","Go to Homepage",function(){d.path("/")}),g.add("g b","Go to Browser",function(){d.path("/browser")}),g.add("g e","Go to Editor",function(){d.path("/editor")}),g.add("g t","Go to Trails",function(){d.path("/trails")})}]);