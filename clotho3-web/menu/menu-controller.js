'use strict';

Application.Primary.controller('MenuCtrl', ['$scope', '$location', 'Collector', 'Clotho', 'PubSub', function($scope, $location, Collector, Clotho, PubSub) {
    Clotho.get('menu_items').then(function(result) {
        $scope.modes = result;
    });

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
        alert("There was an error: " + rejection);
        //verify this works..
        $location.path(previous);
    });

}]);