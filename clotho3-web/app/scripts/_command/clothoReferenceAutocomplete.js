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
 * ATTRIBUTES
 *
 * bindings
 *
 * @attr ngModel {Model=} binding for query. optional (e.g. just use callback)
 * @attr autocompleteTrigger {Number} Keycode to trigger the autocomplete. Autocomplete will be hidden if this attribute is present, until the keycode is typed //todo
 * @attr autocompleteTriggerInclude {Boolean} Whether the trigger in autocompleteTrigger should be included in the text. Defaults to false //todo
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
 * @attr autocompleteOnBackout {Function=}
 * Called when type backspace and query is empty
 * @attr autocompleteOnEnter {Function=}
 * Called when type enter, and not selecting an autocompletion
 *
 * style config
 *
 * @attr autocompletePopupPosition {String=}
 * Position passed to sharablePopupPosition
 * @attr autocompleteWaitTime {Number}
 * Milliseconds before show autocomplete
 * @attr autocompleteDelimiter {Number} //todo - handle dynamic?
 * Keycode for triggering autocomplete selection. Pressing this key will automatically select the first autocompletion, triggering autocompleteOnSelect. This is not the delimiter for trigger the autocomplete. If none is provided, tokens will only be broken when manually selecting a dropdown suggestion. Will not end if string begins with a single or double quote.
 * Note that there may be weird behavior if you allow autocomplete for strings with spaces
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

		const SPACE_KEY = 32,
					ENTER_KEY = 13,
					ESCAPE_KEY = 27;

		//              backspace tab enter   escape  left  up  right down
		var HOT_KEYS = [8,        9,  13,     27,     37,   38, 39,   40];

		//todo - incorporate reference delimiter -- outside this directive?
		HOT_KEYS.push(ClothoReferenceDelimiter.keycode);

		return {
			restrict: 'A',
			scope: {
				query: '=ngModel',
				autocompleteTrigger: '@?',
				forceVisible: '=?',
				autocompletions : '=?',
				autocompleteHasFocus : '=?',
				autocompleteOnSelect: '&?',
				autocompleteOnKeydown : '&?',
				autocompleteOnQuery : '&?',
				autocompleteOnBackout : '&?',
				autocompleteOnEnter : '&?',
				autocompleteDelimiter : '@?',
				autocompletePopupPosition: '@?',
				autocompleteWaitTime : '@?'
			},
			link: function clothoAutocompleteLink(scope, element, attrs, ngModelCtrl) {

				var localHotkeys = angular.copy(HOT_KEYS);

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
					//wrap in apply so propagates to $parent, and will $digest down
					scope.$apply(function () {
						scope.autocompleteHasFocus = false;
					});
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


				scope.setQueryString = function (string) {
					//if string is empty, just get the current contents
					//example would be paste, after timeout just check contents of element
					string = angular.isEmpty(string) ? scope.query : string;

					//reset
					resetQuery();
					resetMatches();

					//todo - handle spaces?
				};

				// select current selected one, otherwise null
				scope.select = function (activeIdx) {

					var selected = activeIdx > -1 ? scope.autocompletions[activeIdx] : null;

					scope.autocompleteOnSelect({
						$item: selected,
						$query : scope.query
					});

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

				/* QUERY CHANGE LISTENER */

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
				//if no query or autcomplete not currently open, blur
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

				/* EVENT LISTENERS */

				scope.autocompleteDelimiter && localHotkeys.push(scope.autocompleteDelimiter);

				//bind keyboard events from HOT_KEYS + delimiter
				element.bind('keydown', function (evt) {

					scope.autocompleteOnKeydown({$event : evt, $keycode : evt.which});

					//typeahead is open and an "interesting" key was pressed
					if (localHotkeys.indexOf(evt.which) === -1) {
						return;
					}

					//reference delimiter
					if (evt.which === ClothoReferenceDelimiter.keycode) {
						//todo

					}
					// token delimiter
					else if (evt.which === scope.autocompleteDelimiter) {
						if (scope.query == '') {
							//don't allow selection of empty, keep showing placeholder
							//allow default to be prevented
						}
						//if first letter is quote, don't end the token
						else if ( ! checkInQuote(scope.query) ) {
							scope.$apply(function () {
								//if there is one result, select it otherwise null (token is query)
								scope.select(scope.autocompletions.length == 1 ? 0 : -1);
							});
							//space is prevented (by preventDefault below) and placeholder shown again
						}
					}
					//backspace
					else if (evt.which === 8) {
						if (scope.query.length) {
							// return to handle deleting single letter or highlighted text
							// (don't prevent default)
							return;
						} else {
							scope.$apply(function () {
								scope.autocompleteOnBackout();
							});
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

					}
					//right
					else if (evt.which === 39) {

					}
					//enter + tab
					else if (evt.which === 13 || evt.which === 9) {
						//if highlighted dropdown select() it, otherwise we'll run enter callback
						if (scope.activeIdx >= 0) {
							scope.$apply(function () {
								scope.select(scope.activeIdx);
							});
						} else {
							//select() the fragment
							if (scope.query.length) {
								scope.$apply(function () {
									scope.select();
								});
							}
							//submit
							if (evt.which == 13) {
								scope.$apply(function () {
									scope.autocompleteOnEnter()
								});
							}
						}
						scope.$digest();
					}
					//escape
					else if (evt.which === 27) {
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

				//on pasting text, pass to setQuery method
				element.on('paste', function (evt) {
					//copied text only available on clipboard, but inconsistent use and access so just do simple workaround and $timeout then process element
					//don't want to preventDefault() if we're getting it next event loop
					$timeout(function () {
						scope.setQueryString(null, true);
					});
				});

				scope.$on('$locationChangeSuccess', function () {
					//timeout because triggers $digest() - schedule outside angular event loop
					setTimeout(resetActive);
				});

				function clothoAutocompleteBlurHandler (event) {

					//first, check if (1) have focus (2) autocomplete inactive
					if ( scope.autocompleteHasFocus && scope.activeIdx < 0 ) {
						// then check if (1) no event passed OR
						// (2A) not active element AND (2B) didn't click child
						// note - popups will not be children of autocomplete, so don't want to hide everything in that event. can use escape to hide stuff.
						if (angular.isUndefined(event) || ($document[0].activeElement != element[0] && !element[0].contains(event.target)) ) {
							//timeout so can prevent default somewhere else
							$timeout(resetActive);
						}
					}
				}

				//we could bind/rebind on focus events...
				$document.bind('click', clothoAutocompleteBlurHandler);
				scope.$on('$destroy', function() {
					$document.unbind('click', clothoAutocompleteBlurHandler);
				});

				/* INIT */

				//whether the input has focus. need our own check because several child elements should not remove perceived focus, but will unfocus the element when interacted with.
				scope.autocompleteHasFocus = scope.autocompleteHasFocus === true;

				//query needs to have value, whether passed in or empty for init
				scope.query = scope.query || '';

				//set up autocompletions array
				resetMatches();
				element.after($compile(listingEl)(scope));
			}
		};
	});
