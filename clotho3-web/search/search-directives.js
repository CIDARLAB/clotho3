'use strict';

Application.Search.directive('clothoSearchbar', ['Clotho', 'Searchbar', '$location', '$window', function(Clotho, Searchbar, $location, $window) {

    return {
        restrict: 'A',
        replace: true,
        templateUrl: "search/searchbar-container.html",
        controller: function($scope, $element, $attrs) {

            $scope.options = Searchbar.options;
            $scope.log = Searchbar.log;
            $scope.autocomplete = Searchbar.autocomplete;
            $scope.display = Searchbar.display;

            $scope.setQuery = Searchbar.setQuery;
            $scope.submit = Searchbar.submit;
            $scope.execute = Searchbar.execute;

            //functions
            $scope.$watch('display.query', function(newValue, oldValue) {
                $scope.display.autocomplete = !!newValue;
                if (!!newValue) {
                    Clotho.autocomplete($scope.display.query);
                }
            });

            $scope.fullPageLog = function() {
                $location.path("terminal");
                $scope.display.hide('log')
            };

            $scope.pathIsTerminal = function() {
                var regexp = /^\/terminal.*$/;
                return regexp.test($location.path());
            };


            /*** help icons ***/

            $scope.newPage = function() {
                $window.open($window.location.origin, "_blank");
            };

            $scope.newWorkspace = function() {
                $window.open($window.location.origin, "_blank");
            };

            $scope.showMeHow = function() {
                console.log("tutorial");
            };

            $scope.aboutClotho = function() {
                console.log("about clotho");
            };

            $scope.toggleTooltips = function() {
                console.log("tooltips");
            };

            //testing

            $scope.sayTest = function() {
                Clotho.say('This is a test message');
            }


        },
        link: function (scope, element, attrs, controller) {

        }
    }
}]);

Application.Search.directive('clothoSearchbarHelppane', ['Clotho', function(Clotho) {

    return {
        restrict: 'A',
        replace: true,
        template: "",
        controller: function($scope, $element, $attrs) {

        },
        link: function($scope, $element, $attrs) {

        }
    }
}]);

Application.Search.directive('clothoSearchbarAutocomplete', ['Clotho', 'Searchbar', function(Clotho, Searchbar) {

    return {
        restrict: 'A',
        replace: true,
        templateUrl: 'search/autocomplete-partial.html',
        controller: function($scope, $element, $attrs) {
            //inherited, not necessary
            //$scope.display = Searchbar.display;
            //$scope.autocomplete = Searchbar.autocomplete;
        },
        link: function($scope, $element, $attrs) {

        }
    }
}]);

Application.Search.directive('clothoSearchbarLog', ['Clotho', 'Searchbar', '$timeout', function(Clotho, Searchbar, $timeout) {

    return {
        restrict: 'A',
        replace: true,
        templateUrl : "search/log-partial.html",
        controller: function($scope, $element, $attrs) {
            //inherited, not necessary
            //$scope.options = Searchbar.options;
            //$scope.log = Searchbar.log;
        },
        link: function($scope, $element, $attrs) {

        }
    }
}]);

Application.Search.controller('TerminalCtrl', ['$scope', 'Clotho', 'Searchbar', '$location', function($scope, Clotho, Searchbar, $location) {
    $scope.log = Searchbar.log;


}]);
