angular.module('clotho.tokenizer')
/**
 * @name clotho-token
 *
 * @description
 * Given a name and UUID, renders a token which can display more information upon interaction via sharable-popup
 *
 * Can either pass token directly, or tokenModel (will overwrite if not empty), or tokenId (will call Clotho.get() and will overwrite if not empty)
 *
 * //todo - better active handling + $scope propagation and tying to collection
 */
	.directive('clothoToken', function ($parse, $document, Clotho, clothoTokenFactory, ClothoSchemas) {
		return {
			restrict: 'E',
			templateUrl: "views/_command/token.html",
			scope: {
				token : '=?',
        tokenId: '@?',
				tokenModel : '=?',
				tokenName : '@?',
				tokenActivePass : '@?',
				popupPlacement : '@?',
				popupTrigger : '@?',
				onClick : '&?',
				onRemove : '&?'
			},
			link: function clothoTokenLink(scope, element, attrs, ngModelCtrl) {

        scope.$watch('token.model', function (newval) {
          var dirtyType = ClothoSchemas.dirtyDetermineType(newval);
          scope.labelClass = 'label-' + ClothoSchemas.typeToColorClass(dirtyType);
          scope.iconClass = ClothoSchemas.determineSharableIcon(dirtyType);

          if (angular.isDefined(scope.token) && scope.token.isSharable()) {
            scope.labelClass += ' isSharable';
          }
        });

				scope.$watch('tokenModel', function (newval) {
					if (!angular.isEmpty(newval)) {
						scope.token = new clothoTokenFactory(newval);
					}
				});

        scope.$watch('tokenId', function (newval) {
          if (!angular.isEmpty(newval)) {
            scope.tokenName = newval;
            Clotho.get(newval, {mute: true}).then(function (result) {
              if (!angular.isEmpty(result)) {
                scope.token = new clothoTokenFactory(result);
              }
            });
          }
        });

        scope.$watch('tokenActivePass', function (newval) {
          if (angular.isDefined(newval)) {
            scope.token.active(scope.$eval(newval));
          }
        });

				element.on('click', function (evt) {
					scope.onClick({$event : evt, $token : scope.token});
					scope.$apply(function () {
						scope.token.active(!scope.token.active());
					});
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
          element.toggleClass('active', !!newval);
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
