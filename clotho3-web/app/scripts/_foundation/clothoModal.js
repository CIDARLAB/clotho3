angular.module('clotho.clothoDirectives')
/**
 * @name clotho-show
 *
 * @description simple modal service which does not rely on $modal from angular.ui, and allows transclusion of child content, or clotho-show if pass widget ID.
 *
 * @usage
 *
 * You can define arbitrary HTML as content or transclude, a template, or a clotho ID to show. All will be shown if defined.
 *
 * (arbitrary HTML, on scope - not angular compiled)
 *
 * <clotho-modal title="Clotho Help" open="showHelp" content="modalContent"></clotho-modal>
 *
 * (arbitrary HTML, transcluded - angular compiled)
 *
<clotho-modal title="arbitrary content" open="true">
 <p>Any content you want here</p>
 <p>bindings are to {{parentModels}}</p>
</clotho-modal>
 *
 * (Clotho.show(<id>)
 *
 <clotho-modal title="clotho.show" id="<WIDGET_ID>"></clotho-modal>
 *
 * (template - angular compiled)
 *
 * <clotho-modal title="Another Template" template-url="'myUrl.html'"></clotho-modal>
 *
 * @description
 * This implementation is a bit different than UI Bootstrap and Angular Strap in that most of the logic is placed in the modal, rather than the service. Usage is similar, though less programmatic manipulation, in favor of ease of creation and ability to support transclusion.
 */
	.directive('clothoModal', function ($parse, $timeout, $http, $compile, $sce, hotkeys) {

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
				title : '@?',
				content : '=?',
				templateUrl : '=?',
				actions : '=?'
			},
			controller : function ($scope, $element, $attrs) {
				//define on controller so available outside template
				$scope.$close = function (event) {
					if ($scope.open) {
						$scope.open = false;
						hotkeys.del('esc');

						$timeout(function () {
							$scope.onClose();
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

				//todo - this should be compiled
				scope.$watch('content', function(newValue, oldValue) {
					scope.contentTrusted = $sce.trustAsHtml(newValue);
				});

				scope.$watch('templateUrl', function (newval, oldval) {
					if (!!newval && (newval != oldval || !oldval || !scope.hasTemplate)) {
						$http.get(newval, {cache: true}).success(function (data, headers) {
							angular.element(element[0].querySelector('[template-insert]')).html($compile(data)(scope));
							scope.hasTemplate = true;
						})
						.error(function (data) {
							//let's hide the area if the template doesn't exist
							scope.hasTemplate = false;
						});
					}
				});

				scope.$watch('open', function (newval, oldval) {
					if (!!newval) {
						scope.open = true;
						hotkeys.add('esc', scope.$close);
						angular.isFunction(scope.onOpen) && scope.onOpen();
					}
				});
			}
		}
	})

/**
 * @name $clothoModal
 *
 * @description
 * Simple service for programmatic creation of clotho-modals.
 *
 * It is likely you may want to use the $modal service provided by angularUI-bootstrap, this is primarily for the minimal builds (API and command bar only), though completely usable elsewhere.
 *
 * This service allows on modal at a time.
 *
 * To pass functions etc, you must create them on the scope, and pass in the $scope. The config options are simply mapped to the DOM as attrs, so should match docs for clotho-modal directive. They should be snake-case. Note that the $scope of the directive is isolate, so all actions must be within the 'scope' of the $scope you pass in.
 *
 * @example (simple)
 * given var msg = <string>, because the string is inserted into the DOM
 *
 *  $clothoModal.create({
			title : "Clotho Alert",
			content : "msg"
		});
 *
 * @example (further demonstrate quoting)
 *
 *  $clothoModal.create({
			title : 'Clotho Login',
			'template-url' : "'views/_command/simpleLogin.html'"
		});
 *
 * @example (passing in scope)
		var newScope = $scope.$new();
		newScope.modalActions =  [
			{
			class : 'info',
			text : 'Great!',
			action : $clothoModal.destroy
			}
		];

	  $clothoModal.create({
			title : "Hey there",
			content : "'here is some <br>content'",
			actions : 'modalActions'
		}, newScope);


 */
	.service('$clothoModal', function($window, $rootScope, $compile) {

		var bodyElement = angular.element($window.document.body);
		var extantModal = null, extantScope = null;

		/**
		 * @name $clothoModal.create
		 *
		 * @description
		 * Create a new clotho-modal element on the body, removing the last one if it existed.
		 *
		 * @param config {Object} directly mapped to attrs of clotho-modal element
		 * @param scope {Scope=} for evaluating scope of modal, otherwise new scope created
		 */
		this.create = function (config, scope) {
			destroy();

			var options = angular.extend({}, config);
			extantScope = scope = scope || $rootScope.$new();

			//reset the close
			var oldClose = scope.$eval(options.onClose) || angular.noop;
			scope.clothoModalClose = function () {
				oldClose();
				destroy();
			};
			options['on-close'] = 'clothoModalClose()';

			extantModal = $compile(angular.element('<clotho-modal>').attr(options))(scope);
			bodyElement.append(extantModal);
		};

		/**
		 * @name $clothoModal.destroy()
		 *
		 * @description
		 * Remove a modal created by this service, otherwise does nothing.
		 */
		function destroy () {
			extantScope && extantScope.$destroy();
			extantModal && extantModal.remove();
			extantModal = null;
			extantScope = null;
		}

		this.destroy = destroy;

	});