/**
 field types that are handled:
 backend: CSS (url), mixin (array|url), script (array|url), onload (array|url)
 content: text (string|html), video (object), template (url), quiz (object), markdown (text), wiki (text)

 you can add a controller using the ng-controller directive in a template, and declare it as a mixin.
 */
angular.module('clotho.trails').directive('trailPage', function ($timeout, $q, $controller, hotkeys) {

	return {
		restrict: 'A',
		templateUrl: 'views/_trails/trailPage.html',
		scope: {
			page: '=trailPage',
			next: '=',
			prev: '='
		},
		link: function trailPageLink(scope, element, attrs) {

			scope.helpModalOpen = false;
			function toggleHelpModal() {
				scope.helpModalOpen = !scope.helpModalOpen;
			}

			scope.$watch('page', function (newPage) {
				if (!angular.isEmpty(newPage)) {
					scope.createPage()
				}
			});

			//future - this provides the foundation for Clotho.view() -- move it there
			scope.createPage = function () {
				if (angular.isDefined(scope.page.dictionary)) {
					angular.extend(scope, scope.page.dictionary);
				}

				return $q.all([
					$clotho.extensions.css(scope.page.css),
					$clotho.extensions.mixin(scope.page.mixin),
					$clotho.extensions.script(scope.page.script)
				])
				.then(function () {
					scope.pageComponents = scope.page.contents;

					if (scope.page.help) {
						hotkeys.add('h', 'Toggle Page Help', toggleHelpModal, false);
					}

					$timeout(function () {
						return $clotho.extensions.script(scope.page.onload)
					});
				});
			};
		}
	}
});