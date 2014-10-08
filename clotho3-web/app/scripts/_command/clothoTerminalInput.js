angular.module('clotho.commandbar')
	//todo - allow delegation of options in clotho-reference-input
	.directive('clothoTerminalInput', function (Clotho, ClientAPI) {

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
					ClientAPI.say({
						from: 'client',
						channel: 'submit',
						class: 'info',
						text: scope.query
					});
					Clotho.submit(scope.query).then(function (response) {
						ClientAPI.say({
							from: 'server',
							channel: 'submit',
							class: 'success',
							text: response
						});
						scope.query = '';
					});
				};
			}
		};
	});
