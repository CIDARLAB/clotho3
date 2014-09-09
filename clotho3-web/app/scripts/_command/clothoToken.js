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
 * //todo - probably just want to pass in collection?
 * //todo - beter active handling + $scope propagation
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
					scope.token.active(scope.$eval(newval));
				});

				scope.$watch('tokenModel', function (newval) {
					if (!angular.isEmpty(newval)) {
						scope.token = new clothoTokenFactory(newval);
					}
				});

				element.on('click', function (evt) {
					scope.$apply(function () {
						scope.token.active(!scope.token.active());
					});
					scope.onClick({$event : evt, $token : scope.token});
				});

        var tokenKeypressListener = function (evt) {

          //backspace
          if (evt.which === 8) {
            evt.preventDefault();
            element.remove();
            scope.$destroy();
          }
          //escape
          else if (evt.which === 27) {
            evt.preventDefault();
            scope.$apply(function() {
              scope.token.active(false);
            });
          }
        };

        scope.$watch('token.active()', function (newval) {
          if (newval) {
            $document.on('keydown', tokenKeypressListener);
          } else {
            $document.off('keydown', tokenKeypressListener);
          }
        });

        scope.$on('$destroy', function (evt) {
          scope.onRemove({$token : scope.token});
          $document.off('keydown', tokenKeypressListener);
        });
			}
		}
	});
