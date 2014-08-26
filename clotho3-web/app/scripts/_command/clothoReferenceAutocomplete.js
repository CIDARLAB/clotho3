angular.module('clotho.tokenizer')
/**
 * @ngdoc directive
 * @name clothoReferenceAutocomplete
 *
 * @description
 * As a directive on an input element, renders an autocomplete with typeahead
 *
 * uses directives clothoAutocompleteListing to create the listing
 *
 * bindings
 *
 * @attr ngModel {Model=} binding for query. optional (e.g. just use callback)
 * @attr forceVisible {Boolean=} if true, force open. if false, force hidden.
 * @attr autocompletions {Array=} Bind to the list of autocompletions
 * @attr autocompleteHasFocus {Boolean=} Bind to whether the input has focus
 *
 * event
 *
 * @attr autocompleteOnSelect {Function=}
 * passed $item, $query (current query string)
 * @attr autocompleteOnQuery {Function=}
 * when query changes. passed $query (new query string) and $old (old query string)
 * @attr autocompleteOnKeydown {Function=}
 * passed $event and $keycode
 *
 * style config
 *
 * @attr autocompletePopupPosition {String=}
 * Position passed to sharablePopupPosition
 * @attr autocompleteWaitTime {Number}
 * Milliseconds before show autocomplete
 * @attr autocompleteTrigger {Number} Keycode for triggering autocomplete. Defaults to reference delimiter in ClothoReferenceDelimiter //todo
 * @attr autocompleteAllowSpaces {Boolean=} Allow spaces in query. This may lead to weird autocompleting behavior. By default, a space will select the current token (or first if none selected). //todo
 *
 * @example //todo
 <input type="text"
 clotho-autocomplete
 ng-model="query"
 placeholder="{{placeholder}}"
 autocomplete-on-query="unsetTokenCollectionActive()"
 autocomplete-on-select="addToken($item, $query)">

 It is the responsibility of other directives to deal with token collection
 */
	.directive('clothoReferenceAutocomplete', function ($q, $parse, $timeout, $compile, $filter, $document, Clotho, ClothoReferenceDelimiter) {

		//              backspace tab enter   escape  left  up  right down
		var HOT_KEYS = [8,        9,  13,     27,     37,   38, 39,   40];
		const SPACE_KEY = 32,
					ENTER_KEY = 13,
					ESCAPE_KEY = 27;

		var tokenDelimiterCode = SPACE_KEY;
		var tokenDelimiterValue = ' ';

		//todo - incorporate reference delimiter

		return {
			restrict: 'A',
			scope: {
				query: '=?ngModel',
				autocompleteTrigger: '@?',
				forceVisible: '=?',
				autocompletions : '=?',
				autocompleteHasFocus : '=?',
				autocompleteOnSelect: '&?',
				autocompleteOnKeydown : '&?',
				autocompleteOnQuery : '&?',
				autocompletePopupPosition: '@?',
				autocompleteAllowSpaces : '@?',
				autocompleteWaitTime : '@?'
			},
			link: function clothoAutocompleteLink(scope, element, attrs, ngModelCtrl) {

				if (element[0].nodeName !== 'INPUT') {
					return;
				}

				//add attributes to text area
				element.attr({
					autocomplete : "off",
					autocorrect : "off",
					autocapitalize : "off",
					spellcheck : "false"
				});

				var initialQuoteRegexp = /^['"].*/;

				//pop-up element used to display matches
				var listingEl = angular.element('<clotho-autocomplete-listing></clotho-autocomplete-listing>');
				listingEl.attr({
					autocompletions: 'autocompletions',
					active: 'activeIdx',
					select: 'select(activeIdx)',
					"has-focus": 'autocompleteHasFocus',
					query: 'query',
					"force-visible" : "{{forceVisible}}",
					"passed-placement" : "{{autocompletePopupPosition}}"
				});

				/* set up scope values */
				
				//whether the input has focus
				scope.autocompleteHasFocus = scope.autocompleteHasFocus || false;
				//query needs to have value, whether passed in or empty for init
				scope.query = scope.query || '';


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
					scope.autocompleteHasFocus = false;
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
						if ( !results || !results.length ) {
							resetMatches();
						} else {
							scope.autocompletions = $filter('limitTo')(results, 10);
						}
					});
				};

				//select current one, otherwise select first
				scope.select = function (activeIdx) {

					var selected = activeIdx > -1 ? scope.autocompletions[activeIdx] : scope.autocompletions[0];

					if (selected) {
						scope.autocompleteOnSelect({
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

				//Declare the timeout promise var outside the function scope so that stacked calls can be cancelled later
				var timeoutPromise;

				scope.$watch('query', function (newval, oldval) {

					scope.autocompleteOnQuery({$query: newval, $old : oldval});

					if (!!newval && newval.length) {
						scope.autocompleteHasFocus = true;

						if (scope.autocompleteWaitTime > 0) {
							if (timeoutPromise) {
								$timeout.cancel(timeoutPromise);//cancel previous timeout
							}
							timeoutPromise = $timeout(function () {
								getAutocompletions(newval);
							}, scope.autocompleteWaitTime);
						} else {
							getAutocompletions(newval);
						}
					} else {
						resetMatches();
					}
				});

				//we need to abstract this out so that we can bind/unbind beyond scope of element
				function escapeHandler () {
					if (!scope.query.length || !scope.autocompletions.length) {
						element[0].blur();
						resetActive();
					} else {
						resetMatches();
						element[0].focus();
						scope.$digest();
					}
				}

				//bind keyboard events from HOT_KEYS + delimiter
				element.bind('keydown', function (evt) {

					scope.autocompleteOnKeydown({$event : evt, $keycode : evt.which});

					//keep delimiter out of HOT_KEYS check because space is a weird default hotkey
					//if type space and not in quote, and only 1 result, will choose it (enter will not)

					//todo - update selecting this way

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
							// return to handle deleting single letter or highlighted text
							// (don't prevent default)
							return;
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
						console.log('escape pressed');
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
						scope.autocompleteHasFocus = true;
					});
				});

				element.on('blur', function (event) {
					if (scope.autocompleteHasFocus) {
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
					if ( scope.autocompleteHasFocus && (angular.isUndefined(scope.tokenCollection) || !scope.tokenCollection.isActive()) && scope.activeIdx < 0 ) {
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
