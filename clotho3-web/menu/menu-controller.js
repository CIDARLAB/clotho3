'use strict';

Application.Primary.controller('MenuCtrl', ['$scope', '$location', '$timeout', 'Collector', 'PubSub', function($scope, $location, $timeout, Collector, PubSub ) {

    $scope.modes = {"items": [
        {"name" : "Editor", "path" : "/editor"},
        {"name" : "Trails", "path" : "/trails"},
        {"name" : "Browser", "path" : "/browser"},
        {"name" : "Plasmid", "path" : "/plasmid"},
        {"name" : "Construction", "path" : "/construction"},
        {"name" : "Schemas", "path":"/schemas"},
        {"name" : "Scripts", "path":"/scripts"}
    ]};

    $scope.$watch(function () {
        return $location.path();
    }, function (newValue, oldValue) {
        if (!$scope.modes) return;
        angular.forEach($scope.modes.items, function(mode, num) {
            var regexp = new RegExp('^' + mode.path + '.*$', ['i']);
            if (regexp.test(newValue)) {
                mode.class = "active";
            } else {
                //remove class
                mode.class = ""
            }
        });
    });

    //initial check of url, not picked up in watch
    $timeout(function() {
        var path = $location.path();
        angular.forEach($scope.modes.items, function(mode, num) {
            var regexp = new RegExp('^' + mode.path + '.*$', ['i']);
            mode.class = (regexp.test(path)) ? 'active' : '';
        });
    }, 0);

    //do hrefs or this make more sense?
    $scope.goToPage = function(mode) {
        $location.path($scope.modes[mode]);
        //note: $location.path().replace() will avoid creation of page in history
    };

    $scope.clearStorage = function() {Collector.clearStorage()};
    $scope.logListeners = function() {PubSub.logListeners()};

    //route handling
    $scope.$on("$routeChangeStart", function (event, next, current) {
        $scope.status = {
            "text" : "Loading...",
            "class" : "progress-striped active progress-warning"
        };
    });
    $scope.$on("$routeChangeSuccess", function (event, current, previous) {
        $scope.status = {
            "text" : "Successfully Loaded!",
            "class" : "progress-success"
        }
    });
    $scope.$on("$routeChangeError", function (current, previous, rejection) {
        $scope.status = {
            "text" : "Fail to Load!",
            "class" : "progress-error"
        };
    });

}]);