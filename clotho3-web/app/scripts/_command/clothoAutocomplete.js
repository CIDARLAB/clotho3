angular.module('clotho.tokenizer')
/**
 * @description
 * As a directive on an input element, renders an autocomplete with typeahead
 *
 * uses directives clothoAutocompleteListing to create the listing
 *
 * tokenCollection is not necessary
 *
 * @attr tokenCollection will update token collection selection etc. if passed in
 * @attr autocompleteOnSelect {Function} passed $item, $query
 *
 * @example
 <input type="text"
    clotho-autocomplete
    ng-model="query"
    placeholder="{{placeholder}}"
    token-collection="tokenCollection"
    autocomplete-on-select="addToken($item, $query)">

 todo - do not rely on ngModel being query
 todo - pass in popup position
 todo - isolate scope

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

		return {
			restrict: 'A',
			require : '?ngModel',
			link: function clothoAutocompleteLink(scope, element, attrs, ngModelCtrl) {

				//add attributes to text area
				element.attr({
					autocomplete : "off",
					autocorrect : "off",
					autocapitalize : "off",
					spellcheck : "false"
				});

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
					angular.isDefined(scope.tokenCollection) && scope.tokenCollection.unsetActive();
					scope.$digest();
				}

				function checkInQuote (query) {
					return initialQuoteRegexp.test(query.charAt(0));
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
						angular.isDefined(scope.tokenCollection) && scope.tokenCollection.unsetActive();

						/*
						//note - check if query value last char is space (blocked for user input) and select if is
						//relevant when typing out programmatically
						if (newval.charAt(newval.length - 1) == tokenDelimiterValue && !checkInQuote(newval)) {
							console.log('last char is space');
							scope.select();
							return;
						}
						*/

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

				scope.setQueryString = function (string, noReset) {

					string = (angular.isUndefined(string) || angular.isEmpty(string)) ? scope.query : string;

					//reset
					resetQuery();
					resetMatches();
					if (!noReset) {
						angular.isDefined(scope.tokenCollection) && scope.tokenCollection.removeAll();
					}

					angular.forEach(string.split(tokenDelimiterValue), function (token) {
						//want to call parent's add token so updates completeQuery as well
						token.length && scope.addToken(token);
					});
				};

				scope.select = function (activeIdx) {

					console.log('selecting ', activeIdx);

					var selected = activeIdx > -1 ? scope.autocompletions[activeIdx] : scope.query;

					if (selected) {
						onSelectCallback(scope, {
							$item: selected,
							$query : scope.query
						});
					}

					resetQuery();

					//return focus to the input element if a match was selected via a mouse click event
					//need scope to update whether focused
					// use timeout to avoid $rootScope:inprog error
					$timeout(function() {
						//reset matches in here so dropdown is within element (for click listener)
						resetMatches();
						element[0].focus();
					}, 0, false);
				};

				//we need to abstract this out so that we can bind/unbind beyond scope of element
				function escapeHandler () {
					if (!scope.query.length || !scope.autocompletions.length) {
						element[0].blur();
						resetActive();
					} else {
						resetMatches();
						element[0].focus();
						angular.isDefined(scope.tokenCollection) && scope.tokenCollection.unsetActive();
						scope.$digest();
					}
				}

				//bind keyboard events from HOT_KEYS + delimiter
				element.bind('keydown', function (evt) {

					//keep delimiter out of HOT_KEYS check because space is a weird default hotkey
					//if type space and not in quote, and only 1 result, will choose it (enter will not)
					if (evt.which === tokenDelimiterCode) {
						//if there's no query, or it's just a space, set it to empty and let default be prevented
						if (scope.query == '') {
							evt.preventDefault();
						}
						//if first letter is quote, don't end the token
						else if ( ! checkInQuote(scope.query) ) {
							scope.$apply(function () {
								//if there is one result, select it otherwise null (token is query)
								scope.select(scope.autocompletions.length == 1 ? 0 : -1);
							});
							//space is prevented and placeholder shown again
							evt.preventDefault();
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
							if (angular.isDefined(scope.tokenCollection)) {
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
						if (angular.isDefined(scope.tokenCollection) && scope.tokenCollection.isActive()) {
							scope.tokenCollection.setPrevActive();
							scope.$digest();
						} else {
							return;
						}
					}
					//right
					else if (evt.which === 39) {
						if (angular.isDefined(scope.tokenCollection) && scope.tokenCollection.isActive()) {
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
						console.log('hit enter', scope.activeIdx);
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
						scope.$digest();
					}
					//escape
					else if (evt.which === 27) {
						//if no query or autcomplete not currently open, blur
						escapeHandler();
					}

					//at bottom so can return out and continue normal action
					evt.preventDefault();
				});

				var escapeCheckHandler = function (evt) {
					if (evt.which === 27) {
						escapeHandler();
					}
				};

				//$timeout so runs after document click
				element.on('focus', function (event) {
					$document.off('keydown', escapeCheckHandler);
					$timeout(function () {
						scope.hasFocus = true;
					});
				});

				element.on('blur', function (event) {
					if (scope.hasFocus) {
						$document.on('keydown', escapeCheckHandler)
					}
				});

				//on pasting text, break up into tokens (unless quoted) and reset query
				element.on('paste', function (evt) {
					//copied text only available on clipboard, but inconsistent use and access so just do simple workaround and $timeout then process element

					//don't want to prevent the event if we're getting it next event loop
					//evt.preventDefault();

					$timeout(function () {
						scope.setQueryString(null, true);
					});
				});

				scope.$on('$locationChangeSuccess', function () {
					//timeout because triggers $digest()
					setTimeout(resetActive);
				});

				function clothoAutocompleteBlurHandler (event) {
					//only trigger if (1) have focus (2) tokens inactive (3) autocomplete inactive
					if ( scope.hasFocus && (angular.isUndefined(scope.tokenCollection) || !scope.tokenCollection.isActive()) && scope.activeIdx < 0 ) {
						//timeout so can prevent default somewhere else
						if (angular.isUndefined(event) || !element[0].contains(event.target)) {
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
