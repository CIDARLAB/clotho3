/**
 field types that are handled:
 backend: CSS (url), mixin (array|url), script (array|url), onload (array|url)
 content: text (string|html), video (object), template (url), quiz (object), markdown (text), wiki (text)

 you can add a controller using the ng-controller directive in a template, and declare it as a mixin.

 //todo - remove next and prev from needed attrs
 */
angular.module('clotho.trails').directive('trailPage', function ($timeout, $q, $controller, hotkeys, Trails, $clothoModal) {

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
				if (scope.helpModalOpen) {
					$clothoModal.create({
						title : 'Trail Help',
						content : 'page.help'
					});
				} else {
					$clothoModal.destroy();
				}
			}

			scope.$watch('page', function (newPage) {
				if (!angular.isEmpty(newPage)) {
					scope.createPage()
				}
			});

			scope.createPage = function () {
				if (angular.isDefined(scope.page.dictionary)) {
					angular.extend(scope, scope.page.dictionary);
				}

				Trails.downloadDependencies(scope.page.dependencies)
				.then(function (onloadFunction) {
					scope.pageComponents = scope.page.contents;

					if (scope.page.help) {
						hotkeys.add('h', 'Toggle Page Help', toggleHelpModal, false);
					}

					$timeout(function () {
						onloadFunction();
					});
				});
			};
		}
	}
});