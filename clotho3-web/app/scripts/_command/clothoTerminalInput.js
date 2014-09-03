angular.module('clotho.commandbar')
	//todo - allow delegation of options in clotho-reference-input
	.directive('clothoTerminalInput', function (Clotho, ClothoReferenceDelimiter) {

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

				scope.submit = function () {
					Clotho.submit(scope.query);
				}
			}
		};
	});