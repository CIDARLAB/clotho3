'use strict';

Application.Search.directive('clothoSearchbar', ['Clotho', 'Searchbar', '$location', '$window', '$timeout', function(Clotho, Searchbar, $location, $window, $timeout) {

    return {
        restrict: 'A',
        replace: true,
        templateUrl: "/search/searchbar.html",
        controller: function($scope, $element, $attrs) {

            $scope.options = Searchbar.options;
            $scope.log = Searchbar.log;
            $scope.autocomplete = Searchbar.autocomplete;
            $scope.display = Searchbar.display;

            $scope.setQuery = Searchbar.setQuery;
            $scope.submit = Searchbar.submit;
            $scope.execute = Searchbar.execute;

            //functions
            //todo - implement functionality of typeahead directive, but don't rely (don't make angular UI a dependency)
            $scope.$watch('display.query', function(newValue, oldValue) {
                $scope.display.autocomplete = !!newValue;
                if (!!newValue) {
                    Clotho.autocomplete($scope.display.query).then(function(data) {
                        $scope.autocomplete.autocompletions = data;
                    });
                }
            });

            /**** click-outside watcher ***/

            //todo - namespace clickOutside
            $scope.$watch('display.autocomplete', function(newValue, oldValue) {
                if (!!newValue) {
                    //console.log('inactivating autocomplete clickOutside');
                    $scope.$broadcast('clickOutside:$active', $scope.$id)
                } else {
                    //console.log('inactivating autocomplete clickOutside');
                    $scope.$broadcast('clickOutside:$inactive', $scope.$id);
                }
            });

            //todo - fix ugly jQuery hacks
            $scope.currentSelected = 1; //assumes that a.close is present and is first child
            $scope.selectAutoNext = function() {
                if (!$scope.display.autocomplete && $scope.display.query) {
                    $scope.display.show('autocomplete');
                    $scope.currentSelected = 1;
                }


                if ($scope.display.autocomplete && $scope.autocomplete.autocompletions.length) {
                    console.log($scope.currentSelected);
                    $('#clothoSearchbarAutocompleteList li:nth-child('+$scope.currentSelected+')').removeClass('active');

                    if ($scope.currentSelected <= $scope.autocomplete.autocompletions.length)
                        $scope.currentSelected += 1;

                    console.log($scope.currentSelected);


                    var current = $('#clothoSearchbarAutocompleteList li:nth-child('+$scope.currentSelected+')');
                    Searchbar.setQuery(current.scope().item);
                    $scope.display.detail(current.scope().item.uuid);
                    current.addClass('active');
                }
            };
            $scope.selectAutoPrev = function() {
                if ($scope.display.autocomplete && $scope.autocomplete.autocompletions.length) {
                    console.log($scope.currentSelected);

                    $('#clothoSearchbarAutocompleteList li:nth-child('+$scope.currentSelected+')').removeClass('active');
                    if ($scope.currentSelected > 1)
                        $scope.currentSelected -= 1;

                    console.log($scope.currentSelected);

                    var current = $('#clothoSearchbarAutocompleteList li:nth-child('+$scope.currentSelected+')');
                    Searchbar.setQuery(current.scope().item);
                    $scope.display.detail(current.scope().item.uuid);
                    current.addClass('active');
                }
            };

            $scope.fullPageLog = function() {
                $location.path("/terminal");
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

        },
        link: function (scope, element, attrs, controller) {

        }
    }
}]);

Application.Search.controller('TerminalCtrl', ['$scope', 'Clotho', 'Searchbar', '$location', function($scope, Clotho, Searchbar, $location) {
    $scope.log = Searchbar.log;
}]);
