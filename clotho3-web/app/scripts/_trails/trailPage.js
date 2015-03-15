/**
 * @ngdoc directive
 * @name trail-page
 *
 * @description
 * handles inserting a trail page into the page.
 * downloads dependencies, and element which is populated with trail page components
 *
 * you can add a controller using the ng-controller directive in a template, and declare it as a mixin.
 *
 */
angular.module('clotho.trails')
.directive('trailPage', function ($timeout, $q, $controller, hotkeys, Trails, $clothoModal) {

	return {
		restrict: 'A',
		templateUrl: 'views/_trails/trailPage.html',
		scope: {
			page: '=trailPage'
		},
		link: function trailPageLink(scope, element, attrs) {

			scope.helpModalOpen = false;
			scope.closeHelpModal = function () {
				scope.helpModalOpen = false;
			};

			function toggleHelpModal() {
				scope.helpModalOpen = !scope.helpModalOpen;
				if (scope.helpModalOpen && !angular.isEmpty(scope.page.help)) {
					$clothoModal.create({
						title : 'Trail Help',
						content : 'page.help',
						'on-close' : 'closeHelpModal()' //in case close with escape
					}, scope);
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
						angular.isFunction(onloadFunction) && onloadFunction();
					});
				});
			};
		}
	}
});