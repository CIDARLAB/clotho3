angular.module('clotho.commandbar')
/**
 * ClothoTokens are essentially wrappers for clotho sharables, or strings. The expect the fields minimally of name, uuid, schema to be a sharable, or just a string for other keywords
 */

	//todo - pending #216, move from uuid -> id and other fields
	.factory('clothoTokenFactory', function (Clotho) {

		//pass UUID to make object, or just pass value as string
		function ClothoToken (sharable) {
			var self = this;
			self.model = sharable;

			if (this.isSharable()) {
				self.fullSharablePromise = Clotho.get(self.model.uuid, {mute : true}).then(function (data) {
					self.fullSharable = data;
				});
			}
		}

		ClothoToken.prototype.readable = function () {
			//todo - refactor to name pending #216
			return this.model.text || this.model;
		};

		ClothoToken.prototype.id = function () {
			//todo - refactor to name pending #216
			return this.model.uuid || null;
		};

		ClothoToken.prototype.isAmbiguous = function () {
			return angular.isArray(this.model);
		};

		ClothoToken.prototype.isSharable = function () {
			return !this.isAmbiguous() && angular.isDefined(this.model.uuid);
		};

		return ClothoToken;

	})
/**
 * Object to handle a collection of ClothoTokens
 */
	.factory('clothoTokenCollectionFactory', function (clothoTokenFactory) {

		function ClothoTokenCollection (startingTokens) {

			this.tokens = [];
			this.currentTokenIndex = -1;

			//todo - test initial tokens

			if (angular.isArray(startingTokens)) {
				angular.forEach(startingTokens, function (token) {
					this.addToken(token);
				});
			}
		}

		//add a token, pass arguments through to clothoTokenFactory
		ClothoTokenCollection.prototype.addToken = function (sharable) {
			this.tokens.push(new clothoTokenFactory(sharable));
		};

		ClothoTokenCollection.prototype.inRange = function (index) {
			return index > -1 && index < this.tokens.length;
		};

		//get token at given index
		ClothoTokenCollection.prototype.getToken = function (index) {
			return this.tokens[index];
		};

		//return index of token
		ClothoTokenCollection.prototype.indexOf = function (token) {
			return this.tokens.indexOf(token);
		};

		//remove token at given index, return it if removed, otherwise false
		ClothoTokenCollection.prototype.removeToken = function (index) {
			if (this.inRange(index)) {
				return this.tokens.splice(index, 1);
			} else {
				return false;
			}
		};

		//remove all tokens
		ClothoTokenCollection.prototype.removeAll = function () {
			this.tokens.length = 0;
		};

		//remove active token if set and return it, otherwise return false
		ClothoTokenCollection.prototype.removeActiveToken = function () {
			if (this.isActive()) {
				var toReturn =  this.removeToken(this.currentTokenIndex);
				this.unsetActive();
				return toReturn;
			} else {
				return false;
			}
		};

		//set token at given index to be active
		ClothoTokenCollection.prototype.setActive = function (index) {
			if (this.inRange(index)) {
				this.currentTokenIndex = index;
				return index;
			} else {
				return false;
			}
		};

		//set token at last position to be active
		ClothoTokenCollection.prototype.setLastActive = function (index) {
			this.setActive(this.tokens.length - 1);
		};

		//set previous token active, based on current active, otherwise last
		ClothoTokenCollection.prototype.setPrevActive = function (index) {
			this.currentTokenIndex = (this.currentTokenIndex > 0 ? this.currentTokenIndex : this.tokens.length) - 1;
		};

		//set next token active, based on current active, otherwise first
		ClothoTokenCollection.prototype.setNextActive = function (index) {
			this.currentTokenIndex = (this.currentTokenIndex + 1) % this.tokens.length;
		};

		//unset active token
		ClothoTokenCollection.prototype.unsetActive = function (index) {
			this.currentTokenIndex = -1;
		};

		//check if token is active, at index when passed, otherwise if any is active
		ClothoTokenCollection.prototype.isActive = function (index) {
			if (angular.isDefined(index)) {
				return this.currentTokenIndex == index;
			} else {
				return this.currentTokenIndex > -1;
			}
		};

		ClothoTokenCollection.prototype.whichActive = function () {
			return this.currentTokenIndex;
		};

		ClothoTokenCollection.prototype.isLastActive = function () {
			return (this.currentTokenIndex === (this.tokens.length - 1))
		};

		return ClothoTokenCollection;
	})
	.directive('clothoTokenizer', function ($parse, clothoTokenCollectionFactory, Debug) {

		var Debugger = new Debug('clothoTokenizer', '#ee7711');

		return {
			restrict: 'E',
			replace: true,
			require: 'ngModel', //avoid isolate scope so model propagates correctly
			templateUrl: "views/_command/tokenizer.html",
			controller: function clothoTokenizerCtrl($scope, $element, $attrs) {

			},
			link: function clothoTokenizerLink(scope, element, attrs, ngModelCtrl) {

				scope.placeholder = attrs.placeholder;

				var startingTags = $parse(attrs.startingTags)(scope);

				var completeQuery = "";

				scope.tokenCollection = new clothoTokenCollectionFactory(startingTags);

				/* updates + watches */

				function updateModel () {
					Debugger.log('updating model (query, tokens)', completeQuery, scope.tokenCollection.tokens);
					ngModelCtrl.$setViewValue({
						query: completeQuery,
						tokens : scope.tokenCollection.tokens
					});
				}

				scope.$watchCollection('tokenCollection.tokens', function () {
					// if tokens changed, query should reflect it.
					// token is either just text, or a full sharable
					// can assume input element is irrelevant
					completeQuery = '';
					//todo - refactor to name pending #216
					angular.forEach(scope.tokenCollection.tokens, function(token) {
						completeQuery  += token.readable() + ' ';
					});

					//update parent model
					updateModel();
				});

				//todo - handle submit / reset to update completeQuery

				/* functionality */

				scope.addToken = function (item) {
					Debugger.log('TOKENIZER_LINK adding token', item);
					scope.tokenCollection.addToken(item);
				};

				scope.removeToken = function (index, model) {
					Debugger.log('TOKENIZER_LINK removing token', index);
					scope.tokenCollection.removeToken(index);
				};

				scope.tokenActive = function (index) {
					return scope.tokenCollection.isActive(index);
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
	.directive('clothoAutocomplete', function (Clotho, $q, $parse, $timeout, $compile, $filter, $document) {

		//              backspace tab enter   escape  left  up  right down
		var HOT_KEYS = [8,        9,  13,     27,     37,   38, 39,   40];
		var SPACE_KEY = 32;
		var ENTER_KEY = 13;
		var tokenDelimiterCode = SPACE_KEY;
		var tokenDelimiterValue = ' ';
		//todo - add attributes (spellcheck, autocapitalize, etc. if necessary)

		return {
			restrict: 'A',
			//require: 'ngModel',
			controller: function clothoAutocompleteCtrl($scope, $element, $attrs) {},
			link: function clothoAutocompleteLink(scope, element, attrs) {

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

				//time to wait before initiating typeahead request
				var waitTime = 0;

				var resetMatches = function() {
					scope.autocompletions = [];
					scope.activeIdx = -1;
				};

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
				scope.query = undefined;

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
					scope.query = '';

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
							scope.$apply(function () {
								scope.query = scope.query.substring(0, scope.query.length - 1);
							});
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
						scope.activeIdx = (scope.activeIdx + 1) % scope.autocompletions.length;
						scope.$digest();

					}
					//up
					else if (evt.which === 38) {
						scope.activeIdx = (scope.activeIdx ? scope.activeIdx : scope.autocompletions.length) - 1;
						scope.$digest();

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
						scope.$apply(function () {
							scope.select(scope.activeIdx);
						});
					}
					//escape
					else if (evt.which === 27) {
						evt.stopPropagation();

						resetMatches();
						scope.tokenCollection.unsetActive();
						scope.$digest();
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

				//can't use 'blur' because will hide list even when item clicked
				//however, don't want to override element.focus() when focused by clicking somewhere in the tokenizerWrap, which will run after element handler due to way events bubble

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
						scope.query = '';
						resetMatches();
					});
				});
				function clothoAutocompleteBlurHandler (event) {
					if (scope.hasFocus) {
						if (!element[0].contains(event.target)) {
							scope.hasFocus = false;
							resetMatches();
							scope.tokenCollection.unsetActive();
							scope.$digest();
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
	})

	/*
	 * internal directive which displays the actual list of autocompletions
	 */
	.directive('clothoAutocompleteListing', function () {
		return {
			restrict:'EA',
			scope:{
				autocompletions:'=',
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
					scope.select({activeIdx:activeIdx});
				};
			}
		};
	})

	/*
	 * internal directive which displays an autocompletion
	 */
	.directive('clothoAutocompleteMatch', function () {
		return {
			restrict:'EA',
			replace: true,
			scope:{
				index:'=',
				match:'=',
				query:'=',
				active : '='
			},
			templateUrl : 'views/_command/autocompleteMatch.html'
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
	.directive('clothoToken', function (Clotho, clothoTokenFactory) {

		//popover or something when hover via clotho.get()

		return {
			restrict: 'E',
			replace: true,
			templateUrl: "views/_command/token.html",
			scope: {
				tokenCollection : '=',
				tokenIndex : '=',
				tokenActive : '=',
				token : '=ngModel',
				onRemove : '&?'
			},
			controller: function clothoTokenCtrl($scope, $element, $attrs) {

			},
			link: function clothoTokenLink(scope, element, attrs, ngModelCtrl) {

				element.on('click', function (evt) {
					//toggle whether token is active
					scope.tokenCollection[scope.tokenActive ? 'unsetActive' : 'setActive'](scope.tokenIndex)
				});

				scope.removeToken = function (evt) {
					evt.preventDefault();
					scope.onRemove({$token : scope.token, $event : evt});
				};
			}
		}
	});