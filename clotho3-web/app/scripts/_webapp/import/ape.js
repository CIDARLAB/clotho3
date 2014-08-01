'use strict';

angular.module('clotho.webapp')
.controller('ApeImportCtrl', function ($scope, Clotho) {

	$scope.converted = [];

	$scope.processGenbankFile = function (file, content, index) {
		Clotho.run('org.andersonlab.py_convertGB', [content]).then(function (r) {
			console.log(r[0]);
			$scope.converted.push(r[0]);
		});
	};

	$scope.createNucSeq = Clotho.create;
});
