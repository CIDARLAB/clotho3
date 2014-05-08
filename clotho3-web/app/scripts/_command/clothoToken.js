angular.module('clotho.tokenizer')
/**
 * @name clotho-token
 *
 * @description
 * Given a name and UUID, renders a token which can display more information upon interaction via sharable-popup
 */
	.directive('clothoToken', function (Clotho, clothoTokenFactory) {
		return {
			restrict: 'E',
			replace: true,
			templateUrl: "views/_command/token.html",
			scope: {
				tokenCollection : '=',
				tokenIndex : '=',
				tokenActive : '=',
				token : '=ngModel',
				onRemove : '&?'
			},
			link: function clothoTokenLink(scope, element, attrs, ngModelCtrl) {

				element.on('click', function (evt) {
					//toggle whether token is active
					scope.tokenCollection[scope.tokenActive ? 'unsetActive' : 'setActive'](scope.tokenIndex)
				});

				scope.removeToken = function (evt) {
					evt.preventDefault();
					scope.onRemove({$token : scope.token, $event : evt});
				};
			}
		}
	});