angular.module('clotho.webapp')
.controller('TeamCtrl', function ($scope, Clotho) {
	$scope.ClothoTeam = Clotho.query({schema : 'LabPerson'});
});
