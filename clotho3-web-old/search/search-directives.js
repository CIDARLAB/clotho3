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
                /*
                //future - reactivate when autocomplete is working (and will need to show autocomplete -- ng-hide)
                $scope.display.autocomplete = !!newValue;
                if (!!newValue) {
                    Clotho.autocomplete($scope.display.query).then(function(data) {
                        $scope.autocomplete.autocompletions = data;
                    });
                }
                */
            });

            //todo - fix ugly jQuery hacks
            //$scope.currentSelected = 1; //assumes that a.close is present and is first child
            $scope.prevSubmittedIndex = false;
            $scope.selectAutoNext = function($event) {
                //temporary - next submitted command
                $event.preventDefault();

                //todo - more robust (doens't really work)

                $scope.prevSubmittedIndex =
                    (!$scope.prevSubmittedIndex) ?
                        0 :
                        (   ($scope.prevSubmittedIndex < $scope.display.queryHistory.length -1) ?
                            $scope.prevSubmittedIndex + 1 :
                            $scope.display.queryHistory.length - 1
                        );

                console.log($scope.prevSubmittedIndex);


                Searchbar.setQuery($scope.display.queryHistory[$scope.prevSubmittedIndex]);


                /*if (!$scope.display.autocomplete && $scope.display.query) {
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
                }*/
            };
            $scope.selectAutoPrev = function($event) {

                //temporary -- previous submitted command
                $event.preventDefault();

                $scope.prevSubmittedIndex =
                    (!!$scope.prevSubmittedIndex) ?
                        (   ($scope.prevSubmittedIndex > 0) ?
                            $scope.prevSubmittedIndex - 1 :
                            0
                        ) :
                        (   ($scope.display.queryHistory.length) ?
                            $scope.display.queryHistory.length - 1 :
                            0
                        );

                console.log($scope.prevSubmittedIndex);

                Searchbar.setQuery($scope.display.queryHistory[$scope.prevSubmittedIndex]);



                /*if ($scope.display.autocomplete && $scope.autocomplete.autocompletions.length) {
                    console.log($scope.currentSelected);

                    $('#clothoSearchbarAutocompleteList li:nth-child('+$scope.currentSelected+')').removeClass('active');
                    if ($scope.currentSelected > 1)
                        $scope.currentSelected -= 1;

                    console.log($scope.currentSelected);

                    var current = $('#clothoSearchbarAutocompleteList li:nth-child('+$scope.currentSelected+')');
                    Searchbar.setQuery(current.scope().item);
                    $scope.display.detail(current.scope().item.uuid);
                    current.addClass('active');
                }*/
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
                Clotho.query({name: 'Learning Clotho'}).then(function (result) {
                    $location.path('/trails/' + result[0].id);
                });
            };

            $scope.goHome = function() {
                $location.path('/');
            };

            $scope.aboutClotho = function() {
                //$window.open('http://www.clothocad.org/index.php/background/', 'aboutClotho');
                $location.path('/about')
            };

            $scope.teamClotho = function() {
                $location.path('/team');
            };

            $scope.toggleTooltips = function() {
                console.log("tooltips");
            };

        },
        link: function (scope, element, attrs, controller) {

        }
    }
}]);