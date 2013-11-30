angular.module('clotho.commandbar').directive('clothoCommandBar', function(Clotho, CommandBar, $location, $window) {

	//todo - implement functionality of typeahead directive, but don't rely (don't make angular UI a dependency)

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
				/*
				 //future - reactivate when autocomplete is working (and will need to show autocomplete -- ng-hide)
				 $scope.display.autocomplete = !!newValue;
				 if (!!newValue) {
				 Clotho.autocomplete($scope.display.query).then(function(data) {
				 $scope.autocomplete.autocompletions = data;
				 });
				 }
				 */
			});

			//todo - rewrite
			//$scope.currentSelected = 1; //assumes that a.close is present and is first child
			$scope.prevSubmittedIndex = false;
			$scope.selectAutoNext = function($event) {
				//temporary - next submitted command
				$event.preventDefault();
				var queriesLen = $scope.display.queryHistory.length -1;


				$scope.prevSubmittedIndex =
					(!$scope.prevSubmittedIndex) ?
						0 :
						(   ($scope.prevSubmittedIndex < queriesLen ) ?
							$scope.prevSubmittedIndex + 1 :
							queriesLen
							);

				console.log($scope.prevSubmittedIndex);


				CommandBar.setQuery($scope.display.queryHistory[$scope.prevSubmittedIndex]);


				/*if (!$scope.display.autocomplete && $scope.display.query) {
				 $scope.display.show('autocomplete');
				 $scope.currentSelected = 1;
				 }


				 if ($scope.display.autocomplete && $scope.autocomplete.autocompletions.length) {
				 console.log($scope.currentSelected);
				 $('#clothoSearchbarAutocompleteList li:nth-child('+$scope.currentSelected+')').removeClass('active');

				 if ($scope.currentSelected <= $scope.autocomplete.autocompletions.length)
				 $scope.currentSelected += 1;

				 console.log($scope.currentSelected);


				 var current = $('#clothoSearchbarAutocompleteList li:nth-child('+$scope.currentSelected+')');
				 CommandBar.setQuery(current.scope().item);
				 $scope.display.detail(current.scope().item.uuid);
				 current.addClass('active');
				 }*/
			};
			$scope.selectAutoPrev = function($event) {

				//temporary -- previous submitted command
				$event.preventDefault();
				var queriesLen = $scope.display.queryHistory.length -1;

				$scope.prevSubmittedIndex =
					(!!$scope.prevSubmittedIndex) ?
						(   ($scope.prevSubmittedIndex > 0) ?
							queriesLen :
							0
							) :
						(   ($scope.display.queryHistory.length) ?
							queriesLen :
							0
							);

				console.log($scope.prevSubmittedIndex);

				CommandBar.setQuery($scope.display.queryHistory[$scope.prevSubmittedIndex]);



				/*if ($scope.display.autocomplete && $scope.autocomplete.autocompletions.length) {
				 console.log($scope.currentSelected);

				 $('#clothoSearchbarAutocompleteList li:nth-child('+$scope.currentSelected+')').removeClass('active');
				 if ($scope.currentSelected > 1)
				 $scope.currentSelected -= 1;

				 console.log($scope.currentSelected);

				 var current = $('#clothoSearchbarAutocompleteList li:nth-child('+$scope.currentSelected+')');
				 CommandBar.setQuery(current.scope().item);
				 $scope.display.detail(current.scope().item.uuid);
				 current.addClass('active');
				 }*/
			};

			$scope.fullPageLog = function() {
				$location.path("/terminal");
				$scope.display.hide('log')
			};

			$scope.pathIsTerminal = function() {
				var regexp = /^\/terminal.*$/;
				return regexp.test($location.path());
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
				//$window.open('http://www.clothocad.org/index.php/background/', 'aboutClotho');
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