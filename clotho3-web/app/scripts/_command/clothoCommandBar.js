/**
 * @name clotho-command-bar
 *
 * @description
 * Creates the Clotho Command Bar
 *
 * @example
 * <div clotho-command-bar></div>
 */

angular.module('clotho.commandbar')
.directive('clothoCommandBar', function(Clotho, CommandBar, $location, $window, $compile, $clothoModal) {

	return {
		restrict: 'A',
		replace: true,
		scope: {},
		templateUrl: "views/_command/commandbar.html",
		controller: function($scope, $element, $attrs) {

			$scope.options = CommandBar.options;
			$scope.log = CommandBar.log;
			$scope.autocomplete = CommandBar.autocomplete;
			$scope.display = CommandBar.display;

			$scope.setQuery = CommandBar.setQuery;
			$scope.submit = CommandBar.submit;
			$scope.execute = CommandBar.execute;

			var showClothoLoginModal = false;
			$scope.toggleLogin = function (force) {
				showClothoLoginModal = angular.isDefined(force) ? force : !showClothoLoginModal;
				if (showClothoLoginModal) {
					$clothoModal.create({
						title : 'Clotho Login',
						'template-url' : "'views/_command/simpleLogin.html'"
					});
				} else {
					$clothoModal.destroy();
				}
			};
		},
		link : function clothoCommandBarLink(scope, element, attrs, controller) {

			/*** help icons ***/

			scope.showMeHow = function() {
				Clotho.query({name: 'Learning Clotho'})
				.then(function (results) {
					Clotho.startTrail(results[0].id);
				});
			};

			scope.goHome = function() {
				$location.path('/');
			};

			scope.aboutClotho = function() {
				$location.path('/about')
			};

			scope.teamClotho = function() {
				$location.path('/team');
			};
		}
	}
});