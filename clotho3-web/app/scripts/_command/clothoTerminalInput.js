angular.module('clotho.commandbar')
	.directive('clothoTerminalInput', function () {

		return {
			restrict: 'E',
			templateUrl: 'views/_command/terminalInput.html',
			scope: {
				placeholder: '@',
				autocompleteTrigger: '@?'
			},
			link : function (scope, element, attrs) {
				scope.focusInput = function () {
					element.find('[clotho-reference-autocomplete]').focus();
				};

				scope.soFar = '';

				scope.selectAutocompletion = function (item, query) {
					console.log('you selected for ' + query, item);
					scope.soFar += query;
				};

				scope.trimSoFar = function () {
					scope.soFar = scope.soFar.substring(0, scope.soFar.length - 1);
				}
			}
		};
	});