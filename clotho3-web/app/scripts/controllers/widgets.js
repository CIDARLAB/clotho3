'use strict';

angular.module('clotho.webapp')
  .controller('WidgetsCtrl', function ($scope, $compile) {

		$scope.myModel = "WHATS up";

		$scope.bootstrapCallback = function(element) {
			console.log('widget controller callback! passed element:', element);
		};

		$scope.bootstrapNewApp = function () {
			console.log('bootstrap start');
			var newEl = angular.element('<clotho-show id="123456789" callback="bootstrapCallback"></clotho-show>');
			angular.element(document).find('insert-widgets-here').append(newEl);
			$compile(newEl)($scope);
			console.log('done compiling');
		};

		$scope.bootstrapNewAppExternal = function () {
			var newEl = angular.element('<clotho-show id="123456789"></clotho-show>');
			angular.element(document).find('#clothoAppWidgets').append(newEl);
			$compile(newEl)($scope);
		};

		$scope.bootstrapSimple = function () {
			var newEl = angular.element('<clotho-show id="813579135"></clotho-show>');
			angular.element(document).find('insert-widgets-here').append(newEl);
			$compile(newEl)($scope);
		};

		$scope.someValue = 'from parent controller not passed down';
  });
