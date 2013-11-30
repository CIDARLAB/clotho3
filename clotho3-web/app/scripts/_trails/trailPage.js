/**
 field types that are handled:
 backend: CSS (url), mixin (array|url), script (array|url), onload (array|url), controller (name, must be mixed in)
 content: text (string|html), video (object), template (url), quiz (object), markdown (text), wiki (text)
 */
angular.module('clotho.trails').directive('trailPage', function($timeout, $q, $controller) {

	return {
		restrict: 'A',
		template: '<div ng-repeat="comp in pageComponents">' +
			'<div trail-page-component="comp"></div>' +
			'</div>',
		scope: {
			page: '=trailPage'
		},
		compile: function compile(tElement, tAttrs, transclude) {
			return {
				pre: function preLink(scope, element, attrs) {

					console.log(scope.page);

					//future in ng-1.2.x, use notify callbacks for updates
					//future - this provides the foundation for Clotho.view() -- move it there

					scope.setPage = function(thePage) {
						if (!!thePage && angular.isObject(thePage))
							scope.page = thePage;
					};

					scope.$watch('page', function (newPage) {
						if (!!newPage) {
							scope.createPage()
						}
					});

					scope.createPage = function () {
						if (!!scope.page.dictionary) {
							angular.extend(scope, scope.page.dictionary);
						}

						return $clotho.extensions.css(scope.page.css)
						.then(function() {
							return $clotho.extensions.mixin(scope.page.mixin)
						})
						.then(function() {
							return $clotho.extensions.script(scope.page.script)
						})
						.then(function (){


							// verify this is the best way to do this
							// todo - create empty as default and use ng-controller in template? -- avoid attaching to DOM directly as may change
							//check for controller, must be already included (e.g. by mixin)
							if (scope.page.controller) {
								var locals = {};
								locals.$scope = scope;
								var ctrl = $controller(scope.page.controller, locals);
								element.data('$ngControllerController', ctrl);
							}

						})
						.then(function() {
							scope.pageComponents = scope.page.contents;
						}, function(error) {

							//todo - better error handling

							console.log(error);
						})
						.then(function () {
							return (!scope.page.onload) ? $q.when() : $timeout(function() {
								console.log('loading page onload script');
								return $clotho.extensions.script(scope.page.onload)
							});
						});
					};

					//init
					return scope.createPage(scope.page);
				},
				post: function postLink(scope, element, attrs) {
					scope.$watch(attrs.trailPage, function(oldpage, newpage) {
						!!newpage && scope.setPage(newpage);
					}, true)
				}
			}
		}
	}
});