
Application.Search.controller('TerminalCtrl', ['$scope', 'Clotho', 'Searchbar', '$location', function($scope, Clotho, Searchbar, $location) {
    $scope.log = Searchbar.log;
}]);