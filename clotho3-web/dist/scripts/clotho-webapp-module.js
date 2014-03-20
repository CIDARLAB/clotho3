angular.module("clotho.webapp",["clotho.foundation","clotho.interface","clotho.dna","ngSanitize","ngRoute"]),angular.module("clotho.webapp").controller("HomeCtrl",["$scope","$location",function(a,b){a.enterClotho=function(){b.path("/trails/bb99191e810c19729de860fe")},a.enterEugene=function(){b.path("/trails/bb02191e810c19729de860aa")},a.enterRaven=function(){b.path("/trails/bb02191e810c19729de860bb")}}]),angular.module("clotho.webapp").controller("AboutCtrl",["$scope",function(a){a.joinMailingList=function(){alert("this doesnt do anything right now"),a.emailToAdd=null}}]),angular.module("clotho.webapp").controller("TeamCtrl",["$scope","Clotho",function(a,b){a.ClothoTeam=b.query({schema:"LabPerson"})}]),angular.module("clotho.webapp").controller("BrowserCtrl",["$scope","Clotho","$filter",function(a,b,c){b.recent().then(function(b){a.recent_array=b,a.sort(!1)}),a.nav={user:{text:"My Stuff",subtext:"Your Digital Locker",value:"user"},group:{text:"Our Group",subtext:"Your Group's Info",value:"group"},all:{text:"Everything",subtext:"All of Clotho",value:"all"}},a.setCurrent=function(b){a.current=b},a.sort=function(b){if(b){if(a.catSort)return;a.catSort=!0,a.recent=c("categorize")(a.recent_array,"type"),a.recent.Instance=c("categorize")(a.recent.Instance,"schema.name")}else{if(!a.catSort&&angular.isDefined(a.catSort))return;a.catSort=!1,a.recent={entries:angular.copy(a.recent_array)}}},a.base64icon=base64icon}]);var base64icon="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEAAAABACAYAAACqaXHeAAACtElEQVR4Xu2Y3UtqURDFlyYVqQhBQon5koqYiUIElSD951qa+AkGUWDUgxC9KEhqfnfXgOLl3gKP2Xlw5kXO0T2zZ+2Z+eG2NBqNCdbYLCqAVoC2gM6ANZ6B0CGoFFAKKAWUAkqBNVZAMagYVAwqBhWDawwB/TOkGFQMKgYVg4pBxeAaK7A0Bh8fH/H6+orJZIL9/X0Eg0FYLJaZpH8og/v7e+zs7CAWi/313Ve6r8LnV7GWEuDh4QH1eh02m038D4dD+Hw++P1+eaYoqVRK3m9ubuLy8hJWq/XbeluFz+8CGhZgNBohnU7LiSYSCQwGA7y8vMDpdMLj8UjMWq0m72jb29u4uLhAu91GtVrFxsYG4vE43t/fwaS3trYQiURwc3OzsM/5ilu0mw0L0O/3ZbPj8ViS/vj4kBYIBAKyh263i0wmg8PDQ7RaLXQ6HRGKmy2VSmg2m9jd3ZV1/I5Vw/VGfS6a+PT3hgVgb5fLZfHDFmCZ05hEOByWJJl4MplEPp9Hr9ebCTAvHte4XC6cnp5iGZ+/LgBPOJvNSumen5/LKeZyOen1k5MTFItF2O12eL1ePD09SaUcHR3JM+35+Vne087OzqSKlvVpRATDFcAT5wxgb1OA6TMFCIVCqFQq/+yHA/Dq6gqcH9fX17Oq2dvbQzQanfkw4tNI8lxjWABOeFYAT+3g4EA+2ddutxvHx8dgmdOYdKFQkOR40kzu7u4Ob29vcDgcUjmsDg5ArjXq89cFYEAOsNvbWzlRGtuBvcwk542tQQFIAc4FCjKlB4Ug9zlHpjRZ1KcpFJhPkCij8UR/ylbh8397M9wCP5Wo2X5UAL0R0hshvRHSGyGzJ7GZ8ZUCSgGlgFJAKWDmFDY7tlJAKaAUUAooBcyexGbGVwooBZQCSgGlgJlT2OzYSoF1p8AnDSiNnx2jBucAAAAASUVORK5CYII=";angular.module("clotho.webapp").controller("EditorCtrl",["$scope","$routeParams","$location","Clotho",function(a,b,c,d){a.editable=b.id,a.editModePass=!1,a.schemas=[],d.query({schema:"Schema"}).then(function(b){a.schemas=b}),a.allInstances=[],d.query({}).then(function(b){for(var c=[],d=0;d<b.length;d++)"BuiltInSchema"!=b[d].schema&&c.push(b[d]);a.allInstances=c}),a.createNewObject=function(){d.create({schema:a.selected}).then(function(b){console.log(b),a.editable=b,a.editModePass=!0}),a.selected=void 0},a.logSelected=function(){console.log(a.selected)},a.createNewSchema=function(){a.editModePass=!0,a.editable={schema:"ClothoSchema",language:"JSONSCHEMA"}}}]),angular.module("clotho.webapp").controller("TrailsCtrl",["$scope","Clotho",function(a,b){a.trails=[],b.query({schema:"Trail"}).then(function(b){a.trails=b}),a.base64icon=base64icon}]),angular.module("clotho.webapp").controller("MenuCtrl",["$scope","$location","$timeout","Collector","PubSub","$keypress","Clotho","$modal",function(a,b,c,d,e,f,g,h){a.modes=[{name:"Editor",path:"/editor"},{name:"Trails",path:"/trails"},{name:"Browser",path:"/browser"}],a.$watch(function(){return b.path()},function(b){a.modes&&angular.forEach(a.modes.items,function(a){var c=new RegExp("^"+a.path+".*$",["i"]);a.class=c.test(b)?"active":""})}),c(function(){var c=b.path();angular.forEach(a.modes.items,function(a){var b=new RegExp("^"+a.path+".*$",["i"]);a.class=b.test(c)?"active":""})},0),a.goToPage=function(c){b.path(a.modes[c])},a.clearStorage=function(){d.clearStorage()},a.logListeners=function(){e.logListeners()},a.reloadModules=function(){g.emit("reloadModels")},a.loggedIn=!1,a.showLogin=function(){h.login().result.then(function(b){console.log(b),b&&(a.username=b,a.loggedIn=!0)})},a.hideNavBar=!0,a.showMenu=function(){console.log("showing"),a.hideNavBar=!a.hideNavBar},f.on("keydown",{"alt-ctrl-comma":"showMenu()"},a)}]),angular.module("clotho.webapp").controller("TrailCtrl",["$scope","$route","$timeout","Clotho","Trails","$keypress",function(a,b,c,d,e,f){a.id=b.current.params.id,a.trail=b.current.locals.trail,a.activate=function(b){console.log(b),b&&angular.isString(b)&&a.current!=b&&(a.current=b,a.currentPage=e.extractPage(a.trail,a.current))},a.home=function(){a.content=a.trail.description,a.current=void 0},a.favorite=function(){e.favorite(a.id)},a.share=function(){e.share(a.id)},a.next=function(){a.activate(e.calcNextPage(a.trail,a.current))},a.prev=function(){a.activate(e.calcPrevPage(a.trail,a.current))},a.mapIcon=e.mapIcon,f.on("keydown",{"alt-right":"next()","alt-left":"prev()"},a),console.log(b.current.params),a.activate(b.current.params.position||"0-0")}]),angular.module("clotho.webapp").controller("WidgetsCtrl",["$scope","$compile",function(a,b){a.myModel="WHATS up",a.myDNAmodel="aaaagt",a.bootstrapCallback=function(a){console.log("widget controller callback! passed element:",a)},a.bootstrapNewApp=function(){console.log("bootstrap start");var c=angular.element('<clotho-show id="123456789" callback="bootstrapCallback"></clotho-show>');angular.element(document).find("insert-widgets-here").append(c),b(c)(a),console.log("done compiling")},a.bootstrapNewAppExternal=function(){var c=angular.element('<clotho-show id="123456789"></clotho-show>');angular.element(document).find("#clothoAppWidgets").append(c),b(c)(a)},a.bootstrapSimple=function(){var c=angular.element('<clotho-show id="813579135"></clotho-show>');angular.element(document).find("insert-widgets-here").append(c),b(c)(a)},a.someValue="from parent controller not passed down"}]);