angular.module('clotho.webapp')
.controller('SettingsCtrl', function ($scope, Clotho) {
	Clotho.get('clotho.developer.maxbates', {mute : true})
	.then(function (r) {
		$scope.person = r;
	})
});
