'use strict';

angular.module('clotho.webapp')
  .controller('WidgetsCtrl', function ($scope, $compile) {
		$scope.bootstrapNewApp = function () {
			console.log('bootstrap start');
			var newEl = angular.element('<clotho-show id="123456789"></clotho-show>');
			$compile(newEl)($scope);
			console.log('done compiling');
			angular.element(document).find('insert-widgets-here').append(newEl);
		};

		$scope.bootstrapNewAppExternal = function () {
			var newEl = angular.element('<clotho-show id="123456789"></clotho-show>');
			$compile(newEl)($scope);
			angular.element(document).find('#clothoAppWidgets').append(newEl);
		};

		$scope.someValue = 'from parent controller not passed down';
  });
