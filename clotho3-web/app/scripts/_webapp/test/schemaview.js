'use strict';

angular.module('clotho.webapp')
  .controller('TestSchemaviewCtrl', function ($scope, Clotho, $clothoModal) {

		Clotho.query({name : "Maxwell Bates"}).then(function(results) {
			$scope.retrievedId = results.length ? results[0].id : '';
		});

		$scope.createModal = function () {
			var newScope = $scope.$new();
			newScope.modalActions =  [
				{
					class : 'info',
					text : 'Great!',
					action : $clothoModal.destroy
				}
			];

			$clothoModal.create({
				title : "Hey there",
				content : "'here is some <br>content'",
				actions : 'modalActions'
			}, newScope);
		}
  });
