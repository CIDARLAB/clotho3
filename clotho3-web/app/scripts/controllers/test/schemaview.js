'use strict';

angular.module('clotho.webapp')
  .controller('TestSchemaviewCtrl', function ($scope, Clotho) {

		Clotho.query({name : "maxbates"}).then(function(results) {
			$scope.retrievedId = results.length ? results[0].id : '';
		});

  });
