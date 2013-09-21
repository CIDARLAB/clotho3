Application.Extensions.controller('ClothoTeamCtrl', ['$scope', 'Clotho', function($scope, Clotho) {

    $scope.ClothoTeam = Clotho.query({schema : 'LabPerson'});

}]);