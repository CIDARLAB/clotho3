angular.module('clotho.commandbar')
	.factory('clothoTokenFactory', function (Clotho) {

		//pass UUID to make object, or just pass value as string
		function ClothoToken (value, uuid) {
			this.value = value;
			this.uuid = uuid || undefined;
			this.isSharable = angular.isDefined(uuid);

			if (this.isSharable) {
				this.fullSharablePromise = Clotho.get(this.uuid).then(function (data) {
					this.fullSharable = data;
				});
			}
		}

		//todo - check ambiguous, check valid

		return ClothoToken

	})
	.factory('clothoTokenCollectionFactory', function (clothoTokenFactory) {

		//fixme - incorporate clothoTokenFactory

		function ClothoTokenCollection (startingTokens) {

			//todo - check initial tokens valid
			this.tokens = angular.isArray(startingTokens) ? startingTokens : [];
			this.currentSelectedIndex = -1;
		}

		ClothoTokenCollection.prototype.getToken = function (index) {
			return this.tokens[index];
		};

		ClothoTokenCollection.prototype.indexOf = function (token) {
			return this.tokens.indexOf(token);
		};

		ClothoTokenCollection.prototype.addToken = function (token) {
			this.tokens.push(token);
		};

		ClothoTokenCollection.prototype.removeToken = function (index) {
			return this.tokens.splice(index, 1);
		};

		ClothoTokenCollection.prototype.removeActiveToken = function () {
			if (this.isActive()) {
				var toReturn =  this.tokens.splice(this.currentSelectedIndex, 1);
				this.unsetActive();
				return toReturn;
			} else {
				return false;
			}
		};

		ClothoTokenCollection.prototype.removeAll = function () {
			this.tokens.length = 0;
		};

		ClothoTokenCollection.prototype.setActive = function (index) {
			if (index > -1 && index < this.tokens.length) {
				this.currentSelectedIndex = index;
			}
		};

		ClothoTokenCollection.prototype.setLastActive = function (index) {
			this.setActive(this.tokens.length - 1);
		};

		ClothoTokenCollection.prototype.setPrevActive = function (index) {
			this.currentSelectedIndex = (this.currentSelectedIndex > 0 ? this.currentSelectedIndex : this.tokens.length) - 1;
		};

		ClothoTokenCollection.prototype.setNextActive = function (index) {
			this.currentSelectedIndex = (this.currentSelectedIndex + 1) % this.tokens.length;
		};

		ClothoTokenCollection.prototype.unsetActive = function (index) {
			this.currentSelectedIndex = -1;
		};

		ClothoTokenCollection.prototype.isActive = function (index) {
			if (angular.isDefined(index)) {
				return this.currentSelectedIndex == index;
			} else {
				return this.currentSelectedIndex > -1;
			}
		};

		return ClothoTokenCollection

	})
	.directive('clothoTokenizer', function ($parse, clothoTokenCollectionFactory) {

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

				function updateModel () {
					console.log('updating model', scope.tokenCollection.tokens);
					ngModelCtrl.$setViewValue(scope.tokenCollection.tokens);
					console.log(ngModelCtrl);
				}

				scope.$watchCollection('tokenCollection.tokens', function () {
					console.log('COLLECTION CHANGED');
					updateModel();
				});

				scope.addToken = function (item) {
					console.log('TOKENIZER_LINK adding token', item);
					scope.tokenCollection.addToken(item);
					//updateModel();
				};

				scope.removeToken = function (index, model) {
					console.log('TOKENIZER_LINK removing token', index);
					scope.tokenCollection.removeToken(index);
					//updateModel();
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
	.directive('clothoAutocomplete', function (Clotho, $q, $parse, $timeout, $compile, $filter) {

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

					var selected = scope.queryResults[activeIdx] || scope.query;

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

				element.on('blur', function () {
					scope.hasFocus = false;
					scope.$apply(function () {
						scope.tokenCollection.unsetActive();
					});
				});

				element.on('focus', function () {
					scope.hasFocus = true;
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
					console.log(activeIdx);
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
				model : '=ngModel',
				onRemove : '&?'
			},
			controller: function clothoTokenCtrl($scope, $element, $attrs) {

			},
			link: function clothoTokenLink(scope, element, attrs, ngModelCtrl) {

				element.on('click', function (evt) {
					if (scope.tokenActive) {
						console.log('sharable object', scope.fullSharable);
						scope.tokenCollection.unsetActive(scope.tokenIndex);
					} else {
						scope.tokenCollection.setActive(scope.tokenIndex);
					}
				});

				scope.removeToken = function (evt) {
					evt.preventDefault();
					scope.onRemove({$model : scope.model});
				};

				//todo - styling based on whether ambiguous

				//todo - allow selection for deletion by autocomplete directive

			}
		}
	});