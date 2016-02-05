"use strict";
angular.module("clotho.fullPackage", ["clotho.foundation", "clotho.commandbar", "clotho.dna", "clotho.interface", "clotho.trails", "clotho.construction", "clotho.webapp", "ngSanitize", "ngRoute"]), angular.module("clothoRoot", ["clotho.fullPackage"]).config(["$routeProvider", "$locationProvider", function(a, b) {
  b.html5Mode(!1).hashPrefix("!"), a.when("/", {
    templateUrl: "views/trail-splash.html",
    controller: "TrailSplashCtrl",
    title: "Home",
    hotkeys: [
      ["h", "Show Intro Modal", "showHelp = !showHelp"]
    ]
  }).when("/about", {
    templateUrl: "views/about.html",
    title: "About"
  }).when("/apps", {
    templateUrl: "views/apps.html",
    title: "Browse Apps"
  }).when("/team", {
    templateUrl: "views/team.html",
    controller: "TeamCtrl",
    title: "Team"
  }).when("/trail", {
    templateUrl: "views/trail.html",
    controller: "TrailCtrl",
    reloadOnSearch: !1,
    resolve: {
      trail: ["Clotho", "$q", "$http", "$route", "Trails", function(a, b, c, d, e) {
        if (angular.isUndefined(d.current.params.id)) return b.when(null);
        var f = b.defer();
        return a.get(d.current.params.id).then(function(a) {
          e.compile(a).then(function(b) {
            d.current.$$route.title = a.name, f.resolve(b)
          })
        }, function() {
          c.get("models/org.clothocad.trails.LearningClotho.json").then(function(a) {
            e.compile(a.data).then(function(a) {
              f.resolve(a)
            })
          })
        }), f.promise
      }]
    },
    hotkeys: [
      ["alt+left", "Previous page of Trail", "prev()"],
      ["alt+right", "Next page of Trail", "next()"]
    ]
  }).otherwise({
    redirectTo: "/"
  })
}]).run(["$rootScope", "$route", "$window", "interfaceConfig", function(a, b, c, d) {
  a.$on("$routeChangeSuccess", function(a, b) {
    var d = angular.isDefined(b.$$route) ? b.$$route.title : null;
    c.document.title = "Clotho Trails" + (angular.isDefined(d) ? " | " + d : "")
  }), a.interfaceConfig = d
}]);
