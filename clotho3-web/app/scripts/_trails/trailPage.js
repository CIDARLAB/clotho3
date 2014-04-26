/**
 field types that are handled:
 backend: CSS (url), mixin (array|url), script (array|url), onload (array|url), controller (name, must be mixed in)
 content: text (string|html), video (object), template (url), quiz (object), markdown (text), wiki (text)
 */
angular.module('clotho.trails').directive('trailPage', function($timeout, $q, $controller) {

	return {
		restrict: 'A',
		templateUrl: 'views/_trails/trailPage.html',
		scope: {
			page: '=trailPage',
			next : '=',
			prev : '='
		},
		compile: function compile(tElement, tAttrs, transclude) {
			return {
				pre: function preLink(scope, element, attrs) {

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
						if (angular.isDefined(scope.page.dictionary)) {
							angular.extend(scope, scope.page.dictionary);
						}

						return $q.all([
							$clotho.extensions.css(scope.page.css),
							$clotho.extensions.mixin(scope.page.mixin),
							$clotho.extensions.script(scope.page.script)
						])
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
					//return scope.createPage(scope.page);
				},
				post: function postLink(scope, element, attrs) {

					/*
					scope.$watch(attrs.trailPage, function(newpage, oldpage) {
						!!newpage && scope.setPage(newpage);
					}, true)
					*/
				}
			}
		}
	}
});