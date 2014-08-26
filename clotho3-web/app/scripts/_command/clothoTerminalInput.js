angular.module('clotho.commandbar')
	.directive('clothoTerminalInput', function () {

		return {
			restrict: 'E',
			templateUrl: 'views/_command/terminalInput.html',
			scope: {
				placeholder: '@'
			},
			link : function (scope, element, attrs) {

			}
		};
	});