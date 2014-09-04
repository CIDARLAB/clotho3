angular.module('clotho.tokenizer')
/**
 * @name clotho-token
 *
 * @description
 * Given a name and UUID, renders a token which can display more information upon interaction via sharable-popup
 *
 * Can either pass token directly, or tokenModel (will overwrite)
 *
 * //todo - better handling of arrays (ambiguity)
 */
	.directive('clothoToken', function ($parse, $document, Clotho, clothoTokenFactory) {
		return {
			restrict: 'E',
			templateUrl: "views/_command/token.html",
			scope: {
				token : '=?',
				tokenModel : '=?',
				tokenName : '@?',
				tokenActivePass : '@?',
				popupPlacement : '@?',
				popupTrigger : '@?',
				onClick : '&?',
				onRemove : '&?'
			},
			link: function clothoTokenLink(scope, element, attrs, ngModelCtrl) {

				scope.$watch('tokenActivePass', function (newval) {
					scope.tokenActive = scope.$eval(newval);
				});

				scope.$watch('tokenModel', function (newval) {
					if (!angular.isEmpty(newval)) {
						scope.token = new clothoTokenFactory(newval);
					}
				});

				element.on('click', function (evt) {
					scope.$apply(function () {
						scope.tokenActive = !scope.tokenActive;
					});
					scope.onClick({$event : evt, $token : scope.token});
				});

        var backspaceEventListener = function (evt) {
          //backspace
          if (evt.which === 8) {
            evt.preventDefault();
            //todo - remove token from parent collection, not just with listener
            element.remove();
            scope.$destroy();
          }
        };

        scope.$watch('tokenActive', function (newval) {
          if (newval) {
            $document.on('keydown', backspaceEventListener);
          } else {
            $document.off('keydown', backspaceEventListener);
          }
        });

        scope.$on('$destroy', function (evt) {
          scope.tokenActive = false;
          scope.onRemove({$token : scope.token});
          $document.off('keydown', backspaceEventListener);
        });
			}
		}
	});