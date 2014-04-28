angular.module('clotho.commandbar')
/**
 * ClothoTokens are essentially wrappers for clotho sharables, or strings. The expect the fields minimally of name, uuid, schema to be a sharable, or just a string for other keywords
 */
	.factory('clothoTokenFactory', function (Clotho) {

		//pass UUID to make object, or just pass value as string
		function ClothoToken (sharable) {
			var self = this;
			self.model = sharable;

			if (this.isSharable()) {
				self.fullSharablePromise = Clotho.get(self.model.uuid).then(function (data) {
					self.fullSharable = data;
				});
			}
		}

<<<<<<< HEAD
		ClothoToken.prototype.readable = function () {
			//todo - refactor to name pending #216
			return this.model.text || this.model;
=======
		ClothoToken.prototype.isAmbiguous = function () {
			//todo
		};

		ClothoToken.prototype.isValid = function () {
			//todo
>>>>>>> d5f512eff9e5aa2d9563e8b347d49975c772b22d
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
<<<<<<< HEAD
			this.currentTokenIndex = -1;
=======
			this.currentSelectedIndex = -1;
>>>>>>> d5f512eff9e5aa2d9563e8b347d49975c772b22d

			//todo - test initial tokens

			if (angular.isArray(startingTokens)) {
				angular.forEach(startingTokens, function (token) {
					this.addToken(token);
				});
			}
		}

		//add a token, pass arguments through to clothoTokenFactory
<<<<<<< HEAD
		ClothoTokenCollection.prototype.addToken = function (sharable) {
			this.tokens.push(new clothoTokenFactory(sharable));
=======
		ClothoTokenCollection.prototype.addToken = function () {
			this.tokens.push(new clothoTokenFactory(arguments));
>>>>>>> d5f512eff9e5aa2d9563e8b347d49975c772b22d
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
<<<<<<< HEAD
				var toReturn =  this.removeToken(this.currentTokenIndex);
=======
				var toReturn =  this.removeToken(this.currentSelectedIndex);
>>>>>>> d5f512eff9e5aa2d9563e8b347d49975c772b22d
				this.unsetActive();
				return toReturn;
			} else {
				return false;
			}
		};

		//set token at given index to be active
		ClothoTokenCollection.prototype.setActive = function (index) {
			if (this.inRange(index)) {
<<<<<<< HEAD
				this.currentTokenIndex = index;
=======
				this.currentSelectedIndex = index;
>>>>>>> d5f512eff9e5aa2d9563e8b347d49975c772b22d
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

				scope.tokenCollection = new clothoTokenCollectionFactory(startingTags);

				/* updates + watches */

				function updateModel () {
					Debugger.log('updating model', scope.tokenCollection.tokens);
					ngModelCtrl.$setViewValue(scope.tokenCollection.tokens);
					Debugger.log(ngModelCtrl);
				}

				scope.$watchCollection('tokenCollection.tokens', function () {
					Debugger.log('COLLECTION CHANGED');
					updateModel();
				});

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
		//todo - add attributes (spellcheck, autocapitalize, etc. if necessary)

		return {
			restrict: 'A',
			//require: 'ngModel',
			controller: function clothoAutocompleteCtrl($scope, $element, $attrs) {},
			link: function clothoAutocompleteLink(scope, element, attrs) {

				var onSelectCallback = $parse(attrs.autocompleteOnSelect);

				//pop-up element used to display matches
				var listingEl = angular.element('<clotho-autocomplete-listing></clotho-autocomplete-listing>');
				listingEl.attr({
					matches: 'queryResults',
					active: 'activeIdx',
					select: 'select(activeIdx)',
					"has-focus": 'hasFocus',
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
						scope.tokenCollection.unsetActive();

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

					var selected = activeIdx > -1 ? scope.queryResults[activeIdx] : scope.query;

					if (selected) {
						onSelectCallback(scope, {
							$item: selected
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

				//bind keyboard events: arrows up(38) / down(40), enter(13) and tab(9), esc(27)
				element.bind('keydown', function (evt) {

					//typeahead is open and an "interesting" key was pressed
					if (HOT_KEYS.indexOf(evt.which) === -1) {
						return;
					}

					evt.preventDefault();

					//backspace
					if (evt.which === 8) {
						if (scope.query.length) {
							scope.$apply(function () {
								scope.query = scope.query.substring(0, scope.query.length - 1);
							});
						} else {
							if (scope.tokenCollection) {
								if (scope.tokenCollection.isActive()) {
									scope.tokenCollection.removeActiveToken();
								} else {
									scope.tokenCollection.setLastActive();
								}
							}
							scope.$digest();
						}
					}
					//down
					else if (evt.which === 40) {
						scope.activeIdx = (scope.activeIdx + 1) % scope.queryResults.length;
						scope.$digest();

					}
					//up
					else if (evt.which === 38) {
						scope.activeIdx = (scope.activeIdx ? scope.activeIdx : scope.queryResults.length) - 1;
						scope.$digest();

					}
					//left
					else if (evt.which === 37) {
						scope.tokenCollection.setPrevActive();
						scope.$digest();
					}
					//right
					else if (evt.which === 39) {
						scope.tokenCollection.setNextActive();
						scope.$digest();
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
						scope.$digest();
					}
				});

				//$timeout so runs after document click
				element.on('focus', function (event) {
					$timeout(function () {
						scope.hasFocus = true;
					});
				});

				//can't use 'blur' because will hide list even when item clicked
				//however, don't want to override element.focus() when focused by clicking somewhere in the tokenizerWrap, which will run after element handler due to way events bubble
				function clothoAutocompleteBlurHandler (event) {
					if (scope.hasFocus) {
						if (!element[0].contains(event.target)) {
							scope.hasFocus = false;
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

<<<<<<< HEAD
				element.on('click', function (evt) {
					//toggle whether token is active
					scope.tokenCollection[scope.tokenActive ? 'unsetActive' : 'setActive'](scope.tokenIndex)
				});
=======
				if (scope.tokenCollection) {
					element.on('click', function (evt) {
						if (scope.tokenActive) {
							console.log('sharable object', scope.fullSharable);
							scope.tokenCollection.unsetActive(scope.tokenIndex);
						} else {
							scope.tokenCollection.setActive(scope.tokenIndex);
						}
					});
				}
>>>>>>> d5f512eff9e5aa2d9563e8b347d49975c772b22d

				scope.removeToken = function (evt) {
					evt.preventDefault();
					scope.onRemove({$token : scope.token, $event : evt});
				};
			}
		}
	});