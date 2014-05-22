angular.module('clotho.webapp')
.controller('TeamCtrl', function ($scope, Clotho) {
	Clotho.query({ id: { $regex: 'clotho.developer.*', $options: 'i' } }).then(function (data) {
		$scope.ClothoTeam	= data;
	});
});
