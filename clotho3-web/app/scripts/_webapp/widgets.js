'use strict';

angular.module('clotho.webapp')
  .controller('WidgetsCtrl', function ($scope, ClientAPI) {

		$scope.myObj = {
			myProp : 'aacccggttt'
		};
		$scope.myModel = 'WHATS up';
		$scope.myDNAmodel = 'aaaagt';

		$scope.bootstrapCallback = function(element) {
			console.log('widget controller callback! passed element:', element);
		};

    //todo - use Clotho.show() once working

		$scope.bootstrapNewApp = function () {
      ClientAPI.display('123456789', 'insert-widgets-here');
		};

		$scope.bootstrapNewAppExternal = function () {
      ClientAPI.display('123456789');
		};

		$scope.bootstrapSimple = function () {
      ClientAPI.display('813579135', 'insert-widgets-here');
		};

		$scope.someValue = 'from parent controller not passed down';

		$scope.openCallback = function () {
			console.log('modal openeed');
		};

		$scope.closeCallback = function() {
			console.log('modal closed');
		};
  });
