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
.directive('clothoCommandBar', function(Clotho, CommandBar, $location, $window) {

	return {
		restrict: 'A',
		replace: true,
		templateUrl: "views/_command/commandbar.html",
		controller: function($scope, $element, $attrs) {

			$scope.options = CommandBar.options;
			$scope.log = CommandBar.log;
			$scope.autocomplete = CommandBar.autocomplete;
			$scope.display = CommandBar.display;

			$scope.setQuery = CommandBar.setQuery;
			$scope.submit = CommandBar.submit;
			$scope.execute = CommandBar.execute;

			//todo - incorporate login, use clotho-modal service

			$scope.showClothoLoginModal = false;
			$scope.showLogin = function() {
				$scope.showClothoLoginModal = true;
			};
		},
		link : function clothoCommandBarLink(scope, element, attrs, controller) {

			/*** help icons ***/

			scope.showMeHow = function() {
				Clotho.query({name: 'Learning Clotho'}).then(function (result) {
					$location.path('/trails/' + result[0].id);
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

			scope.toggleTooltips = function() {
				console.log("tooltips");
			};
		}
	}
});