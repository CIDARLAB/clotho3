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
	.directive('clothoModal', function ($parse, $timeout, hotkeys) {

		return {
			restrict : 'E',
			replace : true,
			transclude : true,
			templateUrl : 'views/_foundation/clothoModal.html',
			scope: {
				id : '@?',
				open : '=?',
				onClose : '&?',
				onOpen : '&?',
				title : '@?'
			},
			controller : function ($scope, $element, $attrs) {
				//define on controller so available in template
				$scope.$close = function (event) {
					if ($scope.open) {
						$scope.open = false;
						hotkeys.del('esc');

						$timeout(function () {
							angular.isFunction($scope.onClose) && $scope.onClose();
						});
					}
				};
			},
			link: function (scope, element, attrs, nullCtrl, transclude) {

				//if the open attr is not present just show it initially
				if (!attrs.open) {
					$timeout(function () {
						scope.open  = true;
					})
				}

				scope.$watch('open', function (newval, oldval) {
					if (!!newval) {
						scope.open = true;
						hotkeys.add('esc', scope.$close);
						angular.isFunction(scope.onOpen) && scope.onOpen();
					}
				});
			}
		}
	});