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
 * we are trying to use the same element in different contexts:
 * 1) command bar - selections should tokenize into tokens, saved outside the input. This flow is easily accomodated - see command bar implementation.
 *
 * 2) terminal - autocomplete should only popup on trigger, and selection should add to model. delimiter should not be in text. Ideally, only the input would change - but this isn't possible. This requires we only autocomplete the last word (i.e. break up by spaces) - which introduces ambiguity, but allows for the most natural flow.
 *
 * 3) textarea - keep the symbol, likely to have multiple in the document
 *
 * 4) others - several options are provided to allow flexibility, but because there are so many, it will take a knowledgeable coder to re-wrap this.
 *
 * note - use the attribute autoGrow to handle dynamic width (e.g. if have text or tokens and displaying as inline-block)
 *
 * note - parent should have style `position: relative` so that listing shows up properly
 *
 * ATTRIBUTES
 *
 * bindings
 *
 * @attr ngModel {Model=}
 * binding for query. optional (e.g. just use callback). Will tie into normal ng-model behavior on the input field.
 * @attr autocompleteTrigger {Boolean}
 * Keycode to trigger the autocomplete. Autocomplete will be hidden if this attribute is present, until the keycode is typed
 * @attr autocompleteTriggerInclude {Boolean}
 * Whether the trigger in autocompleteTrigger should be included in the text. Set to true to enable.
 * @attr autocompleteClearOnSelect {Boolean}
 * When `true`, query will be reset when there is a selection, and behavior should be captured in selection callback. defaults to false. `true` is the behavior of the command bar.
 * @attr autocompleteClearOnEnter {Boolean}
 * When `true`, query will be reset on enter(), and behavior should be captured in autocompleteOnEnter callback. defaults to false. `true` is the behavior of the command bar.
 * @attr forceVisible {Boolean=}
 * Force the visibility of the autocomplete
 * if true, force open. if false, force hidden.
 * @attr autocompletions {Array=}
 * Bind to the list of autocompletions
 * @attr autocompleteHasFocus {Boolean=}
 * Bind to whether the input has focus
 * @attr autocompleteFilter {Filter}
 * Relayed to angular's filter filter.
 * @attr autocompleteBlockInput {Boolean}
 * Block further input. Use instead of [disabled] attribute, so other bindings ok.
 * //todo - add styling
 *
 * event
 *
 * @attr autocompleteOnSelect {Function=}
 * passed $item, $query (current query string). Only called when selecting an autocompletion. Use autocompleteOnEnter to handle queries not tied to an autocompletion
 * @attr autocompleteOnQuery {Function=}
 * when query changes. passed $query (new query string) and $old (old query string)
 * @attr autocompleteOnKeydown {Function=}
 * passed $event and $keycode
 * @attr autocompleteOnBackout {Function=}
 * Called when type backspace and query is empty
 * @attr autocompleteOnEnter {Function=}
 * passed $event, $query
 * Called when type enter, and not selecting an autocompletion.
 * Use this callback to handle query strings which are not tied to a dropdown, and to handle submissions. See clothoReferenceTokenizer for example.
 *
 * style config
 *
 * @attr autocompletePopupPosition {String=}
 * Position passed to sharablePopupPosition
 * @attr autocompleteWaitTime {Number}
 * Milliseconds before show autocomplete
 * @attr autocompleteDelimiter {Number}
 * Keycode for triggering autocomplete selection. This is not the delimiter for trigger the autocomplete - works in a similar way to hitting enter. Pressing this key will automatically select the first autocompletion, triggering autocompleteOnSelect.  If none is provided, tokens will only be broken when manually selecting a dropdown suggestion. Will not end if string begins with a single or double quote.
 * Note that there may be weird behavior if you allow autocomplete for strings with spaces
 *
 * It is the responsibility of other directives to deal with token collection
 *
 * @example See clothoReferenceTokenizer
 * @example see clothoTerminalInput
 * @example
 <textarea clotho-reference-autocomplete
					 class="form-control"
					 ng-model="myModel"
					 ng-trim="false"
					 placeholder="enter a bunch of text"
					 autocomplete-trigger="true"
					 autocomplete-trigger-include="true"
					 autocomplete-popup-position="topRight"></textarea>


 */
	.directive('clothoReferenceAutocomplete', function ($q, $parse, $timeout, $compile, $filter, $document, Clotho, ClothoReferenceDelimiter) {

		const SPACE_KEY = 32,
					ENTER_KEY = 13,
					ESCAPE_KEY = 27;

		//              backspace tab enter   escape  left  up  right down
		var HOT_KEYS = [8,        9,  13,     27,     37,   38, 39,   40];

		HOT_KEYS.push(ClothoReferenceDelimiter.keycode);

		return {
			restrict: 'A',
			require: '?ngModel',
			scope: {
				query: '=ngModel',
				autocompleteTrigger: '=?',
				autocompleteTriggerInclude: '=?',
				autocompleteClearOnSelect: '=?',
        autocompleteClearOnEnter: '=?',
				forceVisible: '=?',
				autocompletions : '=?',
        autocompleteBlockInput: '=?',
        autocompleteFilter : '=?',
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

				function extractReference (query) {
					if (!!scope.autocompleteTrigger) {

						//finds the text after the last clothoReference symbol, up to a space, quote, or delimited passed on scope
						var r = new RegExp('.*' +
							ClothoReferenceDelimiter.symbol +
							'(.+?)(?=[\\s\'"' +
							(angular.isDefined(scope.autocompleteDelimiter) ? scope.autocompleteDelimiter : '') +
							']|$)', 'ig').exec(query);

						return angular.isEmpty(r) ? '' : r[1];
					} else {
						return query;
					}
				}

				//pop-up element used to display matches
				var listingEl = angular.element('<clotho-autocomplete-listing></clotho-autocomplete-listing>');
				listingEl.attr({
					autocompletions: 'autocompletions',
					active: 'activeIdx',
					select: 'select(activeIdx)',
					"has-focus": 'autocompleteHasFocus',
					query: 'query',
					"force-visible" : "{{forceVisible}}",
					"trigger-hide" : 'triggerHide',
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
					//uncomment to delete the current query on reset
					//resetQuery();
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
							scope.activeIdx = -1;

              if (scope.autocompleteFilter) {
                results = $filter('filter')(results, scope.autocompleteFilter)
              }

							scope.autocompletions = $filter('limitTo')(results, 10);
						}
					});
				};

				function breakdownQuery (query) {
					return {
						$query: query,
						$auto: extractReference(query)
					};
				}

				scope.setQueryString = function (string) {
					//if string is empty, just get the current contents
					//example would be paste, after timeout just check contents of element

					var queryString = angular.isEmpty(string) ? scope.query : string;

					scope.query = queryString;
					//note - for some reason, just calling this isn't updating scope.query properly..
					/*if (ngModelCtrl) {
						ngModelCtrl.$setViewValue(queryString);
					}*/

					resetMatches();
				};

				// select current selected one, otherwise null
				scope.select = function (activeIdx) {

					var selected = activeIdx > -1 ? scope.autocompletions[activeIdx] : null;
					var breakdown = breakdownQuery(scope.query);

					scope.autocompleteOnSelect(angular.extend({
            $item: selected
          }, breakdown));

					/*if (scope.autocompleteAddOnSelect === true) {
						scope.setQueryString(breakdown.$rest +
							(breakdown.$rest.length ? ' ' : '') +
							(selected.id ? selected.id : selected));
					}*/

					if (scope.autocompleteClearOnSelect === true) {
						resetQuery();
					} else {
						var replacement = '' +
							(scope.autocompleteTriggerInclude === true ? '@' : '') +
							(selected.id ? selected.id : selected);
						scope.query = scope.query.replace('@' + breakdown.$auto, replacement);
					}

					//return focus to the input element if a match was selected via a mouse click event
					// use timeout to avoid $rootScope:inprog error
					$timeout(function() {
						resetMatches();
						element[0].focus();
					}, 0, false);
				};

				/* QUERY CHANGE LISTENER */

				//Declare the timeout promise var outside the function scope so that stacked calls can be cancelled later
				var timeoutPromise;

				scope.$watch('query', function (newval, oldval) {

					var breakdown = breakdownQuery(scope.query),
							toAutocomplete = breakdown.$auto;

					scope.autocompleteOnQuery(angular.extend({
						$old : oldval
					}, breakdown));

					if (toAutocomplete.length) {
						scope.autocompleteHasFocus = true;

						if (scope.autocompleteWaitTime > 0) {
							if (timeoutPromise) {
								$timeout.cancel(timeoutPromise);//cancel previous timeout
							}
							timeoutPromise = $timeout(function () {
								getAutocompletions(toAutocomplete);
							}, scope.autocompleteWaitTime);
						} else {
							getAutocompletions(toAutocomplete);
						}
					} else {
						resetMatches();
					}
				});

				//we need to abstract this out so that we can bind/unbind beyond scope of element
				//if no query or autcomplete not currently open, blur
				function escapeHandler (event) {
          if (scope.autocompleteTrigger) {
            scope.triggerHide = true;
          }
					if (!scope.query.length || !scope.autocompletions.length) {
						element[0].blur();
						resetActive();
					} else {
						resetMatches();
						element[0].focus();
						scope.$digest();
					}
				}

        //style handler
        scope.$watch('autocompleteBlockInput', function (newval) {
          element.toggleClass('input-disabled', !!newval);
        });

				/* EVENT LISTENERS */

				scope.autocompleteDelimiter && localHotkeys.push(scope.autocompleteDelimiter);

				scope.$watch('autocompleteTrigger', function (newval) {
					scope.triggerHide = !!newval;
				});

				//bind keyboard events from HOT_KEYS + delimiter
				element.bind('keydown', function (evt) {

					scope.autocompleteOnKeydown({$event : evt, $keycode : evt.which});

					//reference delimiter (@)
					//hack - need to check shift state and use alternate keycode because @ is shift+2
					if (evt.which === 50 && evt.shiftKey === true) {
						if (scope.autocompleteTrigger) {
							scope.triggerHide = false;
						}

						//todo - other ways to triggerHide

						scope.$digest();
            //do not return, prevent default
					}

					//typeahead is open and an "interesting" key was pressed
					if (localHotkeys.indexOf(evt.which) === -1) {
            //if blocking input, don't let 'uninteresting' keys do anything
						if (scope.autocompleteBlockInput) {
              evt.preventDefault();
            }
            return;
					}


					// token delimiter - select autocompletion
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
								scope.autocompleteOnBackout({$event: evt});
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
					//side arrows don't do anything now, so let's not prevent default
					//left
					else if (evt.which === 37) {
						return;
					}
					//right
					else if (evt.which === 39) {
						return;
					}
					//enter + tab
					else if (evt.which === 13 || evt.which === 9) {
						//if highlighted dropdown select() it, otherwise we'll run enter callback
						if (scope.activeIdx >= 0) {
              scope.$apply(function () {
                scope.select(scope.activeIdx);
              });
            }
            else {

							/*
							//select() the fragment
							//note - now that we support multiple words in input, not doing this
							if (scope.query.length) {
								scope.$apply(function () {
									scope.select();
								});
							}
							 */
							//enter
							if (evt.which == 13) {
								scope.$apply(function () {
                  scope.autocompleteOnEnter(angular.extend({
                    $event: evt
                  }, breakdownQuery(scope.query)));
                  if (scope.autocompleteClearOnEnter === true) {
                    resetQuery();
                  }
								});
							}
						}
						scope.$digest();
					}
					//escape
					else if (evt.which === 27) {
						escapeHandler(evt);
					}

					//at bottom so can return out and continue normal action
					evt.preventDefault();
				});

				var escapeCheckHandler = function (evt) {
					if (evt.which === 27) {
						escapeHandler(evt);
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
						scope.setQueryString();
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
