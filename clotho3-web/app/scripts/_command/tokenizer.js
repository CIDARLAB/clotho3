angular.module('clotho.commandbar')
	.directive('clothoTokenizer', function ($parse) {

		return {
			restrict: 'E',
			replace: true,
			require: 'ngModel',
			scope: {
				placeholder: "@",
				startingTags: '=',
				model : '=ngModel'
			},
			templateUrl: "views/_command/tokenizer.html",
			controller: function clothoTokenizerCtrl($scope, $element, $attrs) {

				$scope.tokens = [];

			},
			link: function clothoTokenizerLink(scope, element, attrs, ngModelCtrl) {

				function updateModel () {
					console.log('updating model', scope.tokens);
					ngModelCtrl.$setViewValue(scope.tokens);
					console.log(ngModelCtrl);
				}

				scope.addToken = function (item) {
					scope.tokens.push(item);
					updateModel();
				};

				scope.removeToken = function (index, model) {
					scope.tokens.splice(index, 1);
					updateModel();
				};

				scope.focusInput = function () {
					element[0].querySelector('.clothoAutocomplete').focus();
				};
			}
		}
	})
/**
 * Renders an autocomplete, given a query
 */
	.directive('clothoAutocomplete', function (Clotho, $q, $parse, $timeout, $compile, $filter) {

		//              backspace tab enter   escape  up  down
		var HOT_KEYS = [8,        9,  13,     27,     38, 40];

		//todo - add attributes (spellcheck, autocapitalize, etc. if necessary)

		return {
			restrict: 'A',
			require: 'ngModel',
			scope: {
				query: '=ngModel',
				onSelect: '&autocompleteOnSelect'
			},
			controller: function clothoAutocompleteCtrl($scope, $element, $attrs) {

			},
			link: function clothoAutocompleteLink(scope, element, attrs, ngModelCtrl) {

				//pop-up element used to display matches
				var listingEl = angular.element('<clotho-autocomplete-listing></clotho-autocomplete-listing>');
				listingEl.attr({
					matches: 'queryResults',
					active: 'activeIdx',
					select: 'select(activeIdx)',
					hasFocus: 'hasFocus',
					query: 'query'
				});

				scope.hasFocus = false;

				//time to wait before initiating typeahead request
				var waitTime = 0;

				var resetMatches = function() {
					scope.queryResults = [];
					scope.activeIdx = -1;
				};

				var getAutocompletions = function (inputValue) {
					var locals = {$viewValue: inputValue};

					Clotho.autocomplete(scope.query).then(function (results) {
						if (!results || !results.length) {
							resetMatches();
						} else {
							scope.queryResults = $filter('limitTo')(results, 10);
						}
				  });

					/*scope.queryResults = ['Alabama', 'Alaska', 'Arizona', 'Arkansas', 'California', 'Colorado', 'Connecticut', 'Delaware', 'Florida', 'Georgia', 'Hawaii', 'Idaho', 'Illinois', 'Indiana', 'Iowa', 'Kansas', 'Kentucky', 'Louisiana', 'Maine', 'Maryland', 'Massachusetts', 'Michigan', 'Minnesota', 'Mississippi', 'Missouri', 'Montana', 'Nebraska', 'Nevada', 'New Hampshire', 'New Jersey', 'New Mexico', 'New York', 'North Dakota', 'North Carolina', 'Ohio', 'Oklahoma', 'Oregon', 'Pennsylvania', 'Rhode Island', 'South Carolina', 'South Dakota', 'Tennessee', 'Texas', 'Utah', 'Vermont', 'Virginia', 'Washington', 'West Virginia', 'Wisconsin', 'Wyoming'];*/
				};

				//we need to propagate user's query so we can higlight matches
				scope.query = undefined;

				//Declare the timeout promise var outside the function scope so that stacked calls can be cancelled later
				var timeoutPromise;

				scope.$watch('query', function (newval, oldval) {
					if (!!newval && newval.length) {
						scope.hasFocus = true;

						if (waitTime > 0) {
							if (timeoutPromise) {
								$timeout.cancel(timeoutPromise);//cancel previous timeout
							}
							timeoutPromise = $timeout(function () {
								getAutocompletions(newval);
							}, waitTime);
						} else {
							getAutocompletions(newval);
						}
					} else {
						resetMatches();
					}
				});

				scope.select = function (activeIdx) {

					var selected = scope.queryResults[activeIdx] || scope.query;

					scope.onSelect({
						$item: selected
					});

					resetMatches();

					scope.query = '';

					//return focus to the input element if a match was selected via a mouse click event
					// use timeout to avoid $rootScope:inprog error
					$timeout(function() {
						element[0].focus();
					}, 0, false);
				};

				//bind keyboard events: arrows up(38) / down(40), enter(13) and tab(9), esc(27)
				element.bind('keydown', function (evt) {

					//typeahead is open and an "interesting" key was pressed
					if (HOT_KEYS.indexOf(evt.which) === -1) {
						return;
					}

					evt.preventDefault();

					if (evt.which === 8) {
						if (scope.query.length) {
							scope.$apply(function () {
								scope.query = scope.query.substring(0, scope.query.length - 1);
							});
						} else {
							//select last token to delete
							//todo
						}
					} else if (evt.which === 40) {
						scope.activeIdx = (scope.activeIdx + 1) % scope.queryResults.length;
						scope.$digest();

					} else if (evt.which === 38) {
						scope.activeIdx = (scope.activeIdx ? scope.activeIdx : scope.queryResults.length) - 1;
						scope.$digest();

					} else if (evt.which === 13 || evt.which === 9) {
						scope.$apply(function () {
							scope.select(scope.activeIdx);
						});

					} else if (evt.which === 27) {
						evt.stopPropagation();

						resetMatches();
						scope.$digest();
					}
				});

				element.bind('blur', function () {
					scope.hasFocus = false;
				});

				//init()
				resetMatches();
				element.after($compile(listingEl)(scope));
			}
		}
	})

	/*
	 * directive which displays the actual list of autocompletions
	 */
	.directive('clothoAutocompleteListing', function () {
		return {
			restrict:'EA',
			scope:{
				matches:'=',
				query:'=',
				active:'=',
				hasFocus: '=',
				select:'&'
			},
			replace:true,
			templateUrl:'views/_command/autocompleteListing.html',
			link:function (scope, element, attrs) {

				scope.isOpen = function () {
					return scope.hasFocus && scope.matches.length > 0;
				};

				scope.isActive = function (matchIdx) {
					return scope.active == matchIdx;
				};

				scope.selectActive = function (matchIdx) {
					scope.active = matchIdx;
				};

				scope.selectMatch = function (activeIdx) {
					console.log(activeIdx);
					scope.select({activeIdx:activeIdx});
				};
			}
		};
	})

	.directive('clothoAutocompleteMatch', function () {
		return {
			restrict:'EA',
			replace: true,
			scope:{
				index:'=',
				match:'=',
				query:'='
			},
			templateUrl : 'views/_command/autocompleteMatch.html',
			link:function (scope, element, attrs) {

				//todo - currently relies on ngSanitize... need to handle SCE version

				//todo - use highlighting filter

			}
		};
	})

	/* for autocomplete list, bold text matching query */
	.filter('clothoAutocompleteHighlight', function() {

		function escapeRegexp(queryToEscape) {
			return queryToEscape.replace(/([.?*+^$[\]\\(){}|-])/g, "\\$1");
		}

		return function(matchItem, query) {
			return query ? matchItem.replace(new RegExp(escapeRegexp(query), 'gi'), '<strong>$&</strong>') : matchItem;
		};
	})

/**
 * Given a name and UUID, renders a token which can display more information upon interaction
 */
	.directive('clothoToken', function (Clotho) {

		//popover or something when hover via clotho.get()

		return {
			restrict: 'E',
			replace: true,
			templateUrl: "views/_command/token.html",
			require: 'ngModel',
			scope: {
				model : '=ngModel',
				onRemove : '&'
			},
			controller: function clothoTokenCtrl($scope, $element, $attrs) {

			},
			link: function clothoTokenLink(scope, element, attrs, ngModelCtrl) {

				element.on('click', function (evt) {
					Clotho.get(scope.model.uuid).then(function () {

					});
				});

				scope.removeToken = function () {
					scope.onRemove({model : scope.model});
				};


				//todo - styling based on whether ambiguous

				//todo - allow selection for deletion by autocomplete directive

			}
		}
	});