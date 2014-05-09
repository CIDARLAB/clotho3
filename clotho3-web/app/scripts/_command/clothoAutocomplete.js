angular.module('clotho.tokenizer')
/**
 * @description
 * As a directive on an input element, renders an autocomplete with typeahead
 *
 * @example
 <input type="text"
    clotho-autocomplete
    ng-model="query"
    placeholder="{{placeholder}}"
    token-collection="tokenCollection"
    autocomplete-on-select="addToken($item, $query)">

 ngModel is not necessary
 */
	.directive('clothoAutocomplete', function (Clotho, $q, $parse, $timeout, $compile, $filter, $document) {

		//              backspace tab enter   escape  left  up  right down
		var HOT_KEYS = [8,        9,  13,     27,     37,   38, 39,   40];
		var SPACE_KEY = 32;
		var ENTER_KEY = 13;

		var tokenDelimiterCode = SPACE_KEY;
		var tokenDelimiterValue = ' ';

		//time to wait before initiating typeahead request
		var waitTime = 0;

		//todo - easy method to submit (not requiring passage of scope.query)

		//todo - add attributes (spellcheck, autocapitalize, etc. if necessary)

		//todo - allow tie-in of query to model via ngModel

		return {
			restrict: 'A',
			require : '?ngModel',
			link: function clothoAutocompleteLink(scope, element, attrs, ngModelCtrl) {

				var initialQuoteRegexp = /^['"].*/;

				var onSelectCallback = $parse(attrs.autocompleteOnSelect);

				//pop-up element used to display matches
				var listingEl = angular.element('<clotho-autocomplete-listing></clotho-autocomplete-listing>');
				listingEl.attr({
					autocompletions: 'autocompletions',
					active: 'activeIdx',
					select: 'select(activeIdx)',
					"has-focus": 'hasFocus',
					query: 'query'
				});

				scope.hasFocus = false;

				function resetQuery () {
					scope.query = '';
				}

				function resetMatches() {
					scope.autocompletions = [];
					scope.activeIdx = -1;
				}

				function resetActive () {
					resetQuery();
					resetMatches();
					scope.hasFocus = false;
					scope.tokenCollection.unsetActive();
					scope.$digest();
				}

				// get Clotho.autocompletions and update results
				// checks for intiial quote, will not autocomplete empty
				var getAutocompletions = function (inputValue) {

					//check for initial quote
					if (initialQuoteRegexp.test(inputValue)) {
						inputValue = inputValue.substring(1);
					}

					//don't autocomplete empty strings
					if (inputValue.length === 0) {
						return;
					}

					var locals = {$viewValue: inputValue};

					//todo - pending #248 use API option
					Clotho.autocomplete(inputValue).then(function (results) {
						//it no results, or query now empty
						if ( !results || !results.length || !scope.query.length ) {
							resetMatches();
						} else {
							scope.autocompletions = $filter('limitTo')(results, 10);
						}
					});

					/*scope.autocompletions = ['Alabama', 'Alaska', 'Arizona', 'Arkansas', 'California', 'Colorado', 'Connecticut', 'Delaware', 'Florida', 'Georgia', 'Hawaii', 'Idaho', 'Illinois', 'Indiana', 'Iowa', 'Kansas', 'Kentucky', 'Louisiana', 'Maine', 'Maryland', 'Massachusetts', 'Michigan', 'Minnesota', 'Mississippi', 'Missouri', 'Montana', 'Nebraska', 'Nevada', 'New Hampshire', 'New Jersey', 'New Mexico', 'New York', 'North Dakota', 'North Carolina', 'Ohio', 'Oklahoma', 'Oregon', 'Pennsylvania', 'Rhode Island', 'South Carolina', 'South Dakota', 'Tennessee', 'Texas', 'Utah', 'Vermont', 'Virginia', 'Washington', 'West Virginia', 'Wisconsin', 'Wyoming'];*/
				};

				//we need to propagate user's query so we can higlight matches
				//this string represents what is in the autocompete input element
				scope.query = '';

				//Declare the timeout promise var outside the function scope so that stacked calls can be cancelled later
				var timeoutPromise;

				scope.$watch('query', function (newval, oldval) {
					if (!!newval && newval.length) {
						scope.hasFocus = true;
						scope.tokenCollection.unsetActive();

						//todo - don't cancel previous until launch new
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

					var selected = activeIdx > -1 ? scope.autocompletions[activeIdx] : scope.query;

					if (selected) {
						onSelectCallback(scope, {
							$item: selected,
							$query : scope.query
						});
					}

					resetMatches();
					resetQuery();

					//return focus to the input element if a match was selected via a mouse click event
					// use timeout to avoid $rootScope:inprog error
					$timeout(function() {
						element[0].focus();
					}, 0, false);
				};

				//bind keyboard events from HOT_KEYS + delimiter
				element.bind('keydown', function (evt) {

					//keep delimiter out of HOT_KEYS check because space is a weird default hotkey
					//if type space and not in quote, and only 1 result, will choose it (enter will not)
					if (evt.which === tokenDelimiterCode) {
						//if first letter is quote, don't end the token
						if ( ! (initialQuoteRegexp.test(scope.query.charAt(0))) ) {
							scope.$apply(function () {
								//if there is one result, select it otherwise null (token is query)
								scope.select(scope.autocompletions.length == 1 ? 0 : -1);
							});
							//return so space is not prevented
							return;
						}
					}

					//typeahead is open and an "interesting" key was pressed
					if (HOT_KEYS.indexOf(evt.which) === -1) {
						return;
					}

					//backspace
					if (evt.which === 8) {
						if (scope.query.length) {
							//return to handle deleting single letter or highlighted text
							return;
						} else {
							if (scope.tokenCollection) {
								if (scope.tokenCollection.isActive()) {
									var previousActive = scope.tokenCollection.whichActive();
									scope.tokenCollection.removeActiveToken();
									scope.tokenCollection.setActive(previousActive);
								} else {
									scope.tokenCollection.setLastActive();
								}
							}
							scope.$digest();
						}
					}
					//down
					else if (evt.which === 40) {
						if (scope.autocompletions.length) {
							scope.activeIdx = (scope.activeIdx + 1) % scope.autocompletions.length;
							scope.$digest();
						}
					}
					//up
					else if (evt.which === 38) {
						if (scope.autocompletions.length) {
							scope.activeIdx = (scope.activeIdx ? scope.activeIdx : scope.autocompletions.length) - 1;
							scope.$digest();
						}
					}
					//left
					else if (evt.which === 37) {
						if (scope.tokenCollection.isActive()) {
							scope.tokenCollection.setPrevActive();
							scope.$digest();
						} else {
							return;
						}
					}
					//right
					else if (evt.which === 39) {
						if (scope.tokenCollection.isActive()) {
							if ( scope.tokenCollection.isLastActive() ) {
								scope.tokenCollection.unsetActive();
							} else {
								scope.tokenCollection.setNextActive();
							}
							scope.$digest();
						} else {
							return;
						}
					}
					//enter + tab
					else if (evt.which === 13 || evt.which === 9) {
						//if highlighted dropdown select it, otherwise we'll submit
						if (scope.activeIdx >= 0) {
							scope.$apply(function () {
								scope.select(scope.activeIdx);
							});
						} else {
							//if there's an open token, close it
							if (scope.query.length) {
								scope.$apply(function () {
									scope.select();
								});
							}
							//submit
							scope.$apply(function () {
								scope.submit();
							});
						}
					}
					//escape
					else if (evt.which === 27) {
						resetMatches();
						scope.tokenCollection.unsetActive();
						scope.$digest();

						//if there is no query, blur
						if (scope.query.length) {
							evt.stopPropagation();
						} else {
							scope.hasFocus = false;
							element.blur();
							scope.$digest();
						}
					}

					//at bottom so can return out and continue normal action
					evt.preventDefault();
				});

				//$timeout so runs after document click
				element.on('focus', function (event) {
					$timeout(function () {
						scope.hasFocus = true;
					});
				});

				//on pasting text, break up into tokens (unless quoted) and reset query
				element.on('paste', function (evt) {
					//copied text only available on clipboard, but inconsistent use and access so just do simple workaround and $timeout then process element

					//don't want to prevent the event if we're getting it next event loop
					//evt.preventDefault();

					//get the value, split into tokens and save, reset query
					$timeout(function () {
						angular.forEach(scope.query.split(tokenDelimiterValue), function (token) {
							//want to call parent's add token so updates completeQuery as well
							scope.addToken(token);
						});
						resetQuery();
						resetMatches();
					});
				});

				scope.$on('$locationChangeSuccess', function () {
					//timeout because triggers $digest()
					setTimeout(resetActive);
				});

				function clothoAutocompleteBlurHandler (event) {
					//only trigger if (1) have focus (2) tokens inactive (3) autocomplete inactive
					if ( scope.hasFocus && !scope.tokenCollection.isActive() && scope.activeIdx < 0 ) {
						//timeout so can prevent default somewhere else
						if (!element[0].contains(event.target)) {
							$timeout(resetActive);
						}
					}
				}

				$document.bind('click', clothoAutocompleteBlurHandler);
				scope.$on('$destroy', function() {
					$document.unbind('click', clothoAutocompleteBlurHandler);
				});

				//init()
				resetMatches();
				element.after($compile(listingEl)(scope));
			}
		}
	});
