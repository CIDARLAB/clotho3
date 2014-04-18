angular.module('clotho.clothoDirectives')
/**
 * @name clotho-show
 *
 * @description simple modal service which does not rely on $modal from angular.ui, and allows transclusion of child content, or clotho-show if pass widget ID.
 *
 * @usage
 *
 * (arbitrary HTML)
 *
<clotho-modal title="arbitrary content" open="true">
 <p>Any content you want here</p>
 <p>bindings are to {{parentModels}}</p>
</clotho-modal>
 *
 * (Clotho.show(<id>)
 *
 <clotho-modal title="clotho.show" id="<WIDGET_ID>">
 </clotho-modal>
 */
	.directive('clothoModal', function ($parse, $timeout, $keypress) {

		return {
			restrict : 'E',
			replace : true,
			transclude : true,
			templateUrl : 'views/_foundation/clothoModal.html',
			scope: {
				id : '@?',
				open : '=?',
				onClose : '=?',
				onOpen : '=?',
				title : '@?'
			},
			link: function (scope, element, attrs, nullCtrl, transclude) {

				//if the open attr is present, tie a watch, otherwise just show at compile
				if (attrs.open) {
					scope.$watch('open', function (newval, oldval) {
						if (!!newval) {
							scope.open = true;
							angular.isFunction(scope.onOpen) && scope.onOpen();
						}
					});
				}
				else {
					scope.open = true;
					angular.isFunction(scope.onOpen) && scope.onOpen();
				}


				//todo - allow buttons to be pressed and stuff, especially for a show

				scope.$close = function (event) {
					if (scope.open) {
						scope.open = false;

						$timeout(function () {
							angular.isFunction(scope.onClose) && scope.onClose();
						});
					}
				};

				$keypress.on('keyup', {'esc' : '$close()'}, scope);
			}
		}
	});