angular.module('clotho.commandbar').directive('clothoCommandBar', function(Clotho, CommandBar, $location, $window, hotkeys) {

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

			//functions

			$scope.$watch('display.query', function(newValue, oldValue) {
				 $scope.display.autocomplete = !!newValue;
				 if (!!newValue && angular.isString(newValue)) {
					 Clotho.autocomplete($scope.display.query).then(function(data) {
					  $scope.autocomplete.autocompletions = data;
					 });
				 }
			});

			//$scope.currentSelected = 1; //assumes that a.close is present and is first child
			$scope.activeIndex = false;
			$scope.selectAutoNext = function($event) {
				$event.preventDefault();
				$scope.activeIndex = ($scope.activeIndex + 1) % $scope.display.queryHistory.length;

				CommandBar.setQuery($scope.display.queryHistory[$scope.activeIndex]);

			};
			$scope.selectAutoPrev = function($event) {
				$event.preventDefault();
				$scope.activeIndex = ($scope.activeIndex ? $scope.activeIndex : $scope.display.queryHistory.length) - 1;

				CommandBar.setQuery($scope.display.queryHistory[$scope.activeIndex]);
			};

			$scope.fullPageLog = function() {
				$location.path("/terminal");
				$scope.display.hide('log')
			};

			$scope.hideAutocomplete = function () {
				$scope.display.hide('autocomplete');
				$scope.display.undetail();
			};

			$scope.pathIsTerminal = function() {
				var regexp = /^\/terminal.*$/;
				return regexp.test($location.path());
			};


			//todo - incorporate login, use clotho-modal

			$scope.showClothoLoginModal = false;
			$scope.showLogin = function() {
				$scope.showClothoLoginModal = true;
			};


		},
		link : function clothoCommandBarLink(scope, element, attrs, controller) {

			/*** help icons ***/

			scope.newPage = function() {
				$window.open($window.location.origin, "_blank");
			};

			scope.newWorkspace = function() {
				$window.open($window.location.origin, "_blank");
			};

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